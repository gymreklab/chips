/*
  This file defines a producer-consumer queue for threading 
*/

#ifndef THREAD_QUEUE_H_
#define THREAD_QUEUE_H_

#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

template <typename T>
class TaskQueue{
  public:
    T pop();
    void push(const T& item);

  private:
    std::queue<T> _queue;
    std::mutex _mutex;
};

#endif
