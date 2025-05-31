#ifndef __RABOTY_H__
#define __RABOTY_H__

#include "SDL2/SDL.h"
#include "SDL2/SDL_thread.h"

#include <stdbool.h>

#ifndef NUM_OF_THREADS // including main thread
    #define NUM_OF_THREADS (4)
#endif

#ifndef MAX_NUMBER_OF_JOBS
    #define MAX_NUMBER_OF_JOBS (1024)
#endif

struct job
{
    void (*function)(void *);
    void *arguments;
};

void jobs_init(void);
bool job_submit(struct job job);
void jobs_complete_all_work(void);
void jobs_shutdown(void);

uint32_t jobs_get_proccessor_count(void);

#endif // __RABOTY_H__

#ifdef RABOTY_IMPLEMENTATION

struct thread_info
{
    SDL_threadID id;
    int          index;
};

struct job_queue
{
    struct job jobs[MAX_NUMBER_OF_JOBS];
    uint32_t volatile read_index;
    uint32_t volatile write_index;
};

struct job_system
{
    struct job_queue job_queue;
    SDL_sem         *job_semaphore;

    uint32_t number_of_jobs;
    uint32_t volatile number_of_jobs_complete;

    size_t             number_of_threads;
    struct thread_info info[NUM_OF_THREADS]; // index 0 will be main thread
};

static struct job_system job_state = {0};

// compares the dest value with the compare value
// if (dest == compare)
//      dest = exchange
// return initial value of dest
static inline uint32_t platform_interlocked_compare_exchange(volatile uint32_t *dest, uint32_t exchange, uint32_t compare)
{
    return __sync_val_compare_and_swap(dest, compare, exchange);
}

static inline uint32_t platform_interlocked_increment(volatile uint32_t *addend)
{
    return __sync_fetch_and_add(addend, 1);
}

static inline void _memory_barrier(void)
{
    __sync_synchronize();
}

#if 1
static bool _do_work_from_queue(SDL_threadID worker_thread_id)
{
    // SDL_UNUSED(worker_thread_id);

    struct job_queue *job_queue = &job_state.job_queue;

    uint32_t current_read_index, new_read_index;
    do
    {
        current_read_index = job_queue->read_index;
        new_read_index     = (current_read_index + 1) % MAX_NUMBER_OF_JOBS;

        if (current_read_index == job_queue->write_index)
        {
            return true; // sleep
        }

    } while (platform_interlocked_compare_exchange(&job_queue->read_index,
                                                   new_read_index,
                                                   current_read_index) != current_read_index);

    struct job job = job_queue->jobs[current_read_index];

    if (job.function)
    {
        // LOGE("[ERROR] Worker %lu got job at index %u\n", worker_thread_id, current_read_index);
        job.function(job.arguments);
    }
    else
    {
        LOGE("[ERROR] Worker %lu got NULL job at index %u\n", worker_thread_id, current_read_index);
    }

    platform_interlocked_increment(&job_state.number_of_jobs_complete);
    return false; // don't sleep
}
#else
// static bool _do_work_from_queue(SDL_threadID worker_thread_id)
//{
//     worker_thread_id = worker_thread_id; // unused

//    struct job_queue *job_queue = &job_state.job_queue;

//    // read the current read position
//    const uint32_t current_read_index = job_queue->read_index;
//    const uint32_t new_read_index     = (current_read_index + 1) % MAX_NUMBER_OF_JOBS;

//    if (current_read_index != job_queue->write_index)
//    {
//        const size_t index = platform_interlocked_compare_exchange(&job_queue->read_index,
//                                                                   new_read_index,
//                                                                   current_read_index);
//        if (index == current_read_index)
//        {
//            struct job job = job_queue->jobs[index];

//            if (job.function)
//            {
//                job.function(job.arguments);

//                platform_interlocked_increment(&job_state.number_of_jobs_complete);
//            }
//            else
//            {
//                fprintf(stdout, "[ERROR] Worker %lu got NULL job at index %u\n", worker_thread_id, current_read_index);
//            }

//        } // else, another thread has beat me here and I should not do anything
//    }
//    else
//    {
//        return true; // sleep
//    }
//    return false; // dont sleep the thread
//}
#endif

static int _worker_thread_function(void *data)
{
    UNUSED(data);
    // const struct thread_info *thread_info = (const struct thread_info *)(data);

    // SDL_threadID id = thread_info->id;
    SDL_threadID id = SDL_GetThreadID(NULL);

    while (1)
    {
        if (_do_work_from_queue(id))
        {
            SDL_SemWait(job_state.job_semaphore);
        }
    }

    return 0;
}

uint32_t jobs_get_proccessor_count(void)
{
    SYSTEM_INFO sysInfo;
    GetSystemInfo(&sysInfo);

    return sysInfo.dwNumberOfProcessors;
}

void jobs_init(void)
{
    job_state.number_of_threads = NUM_OF_THREADS;

    // initialize job semaphore with a count of 0, indicating no jobs are pending
    const uint32_t initial_count = 0;
    job_state.job_semaphore      = SDL_CreateSemaphore(initial_count);

    if (job_state.job_semaphore == NULL)
        return;

    memset(job_state.job_queue.jobs, 0, sizeof(struct job) * MAX_NUMBER_OF_JOBS);
    job_state.job_queue.read_index  = 0;
    job_state.job_queue.write_index = 0;

    LOG("Job system starting additional %d threads\n", NUM_OF_THREADS - 1);

    for (int thread_index = 1; thread_index < NUM_OF_THREADS; thread_index++)
    {
        struct thread_info *thread_info = job_state.info + thread_index;

        // TODO : do we need to store any thread info? even the id?
        SDL_Thread *t = SDL_CreateThread(_worker_thread_function, "worker", (void *)(thread_info));
        ASSERT(t, "Thread was not valid\n");

        // SDL_threadID id    = SDL_GetThreadID(t); // NOTE : returns 0 for some reason?
        // thread_info->id    = id;
        thread_info->index = thread_index;

        SDL_DetachThread(t);
    }

    job_state.info[0].id    = SDL_ThreadID();
    job_state.info[0].index = 0;
}

bool job_submit(struct job job)
{
    struct job_queue *queue = &job_state.job_queue;

    // trying to write to the last entry to read
    const uint32_t write_index         = queue->write_index;
    const uint32_t next_entry_to_write = (write_index + 1) % MAX_NUMBER_OF_JOBS;

    struct job *new_job = queue->jobs + write_index;

    *new_job = job;

    _memory_barrier();

    job_state.number_of_jobs++;

    queue->write_index = next_entry_to_write;

    SDL_SemPost(job_state.job_semaphore); // signal the worker threads

    return true;
}

// void jobs_complete_all_work(void)
//{
//     while (!_do_work_from_queue(0))
//         ;

//    job_state.number_of_jobs = job_state.number_of_jobs_complete = 0;

//    job_state.job_queue.read_index = job_state.job_queue.write_index = 0;
//}

void jobs_complete_all_work(void)
{
    uint32_t num_jobs      = job_state.number_of_jobs;
    uint32_t complete_jobs = job_state.number_of_jobs_complete;
    do
    {
        _do_work_from_queue(0);

        num_jobs      = job_state.number_of_jobs;
        complete_jobs = job_state.number_of_jobs_complete;
    } while (num_jobs != complete_jobs);

    job_state.number_of_jobs = job_state.number_of_jobs_complete = 0;
    job_state.job_queue.read_index = job_state.job_queue.write_index = 0;
}

void jobs_shutdown(void)
{
    if (job_state.number_of_jobs != 0)
        jobs_complete_all_work();

    memset(job_state.job_queue.jobs, 0, sizeof(struct job) * MAX_NUMBER_OF_JOBS);

    if (job_state.job_semaphore) SDL_DestroySemaphore(job_state.job_semaphore);
}

#endif // RABOTY_IMPLEMENTATION
