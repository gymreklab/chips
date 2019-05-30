#include "src/multithread.h"

template <typename T>
T TaskQueue <T> ::pop(){
  std::unique_lock<std::mutex> mlock(_mutex);
  if (_queue.empty()){
    throw std::out_of_range("Task queue is empty now. All jobs have been finished");
  }
  auto item = _queue.front();
  _queue.pop();
  mlock.unlock();
  return item;
}


template <typename T>
void TaskQueue <T> ::push(const T& item){
  std::unique_lock<std::mutex> mlock(_mutex);
  _queue.push(std::move(item));
  mlock.unlock();
}

template <typename T>
void TaskVector <T> ::push_back(const T& item){
  std::unique_lock<std::mutex> mlock(_mutex);
  _vec.push_back(std::move(item));
  mlock.unlock();
}

