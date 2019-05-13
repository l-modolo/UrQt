/*
UrQt quality and poly nucleotide trimming tool
Copyright (C) 2013  Laurent Modolo

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef DEF_mThreadDone
#define DEF_mThreadDone

#include <future>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <string>
#include "mThreadCircular.hpp"


using namespace std;

template <typename T>
class mThreadDone
{
  public:
  mThreadDone();
  mThreadDone(int size);
  ~mThreadDone();

  void add(T* x, int number);

  void stop();
  void join();

  string output();

  protected:
  inline void get(T** x);
  inline void push_back(T* x, int number);
  inline T* pop_front();

  inline int can_add(int number);
  inline int can_get();
  inline int true_size();
  int mThreadDone_iteration;
  int mThreadDone_iteration_done;

  thread mThreadDone_thread;
  void run();

  // communication between the thread
  mutex mThreadDone_onebyone;
  mutex mThreadDone_empty; // if the list is empty we wait new jobs
  mutex mThreadDone_full; // if the list is full we wait before new jobs
  mutex mThreadDone_initialized;
  condition_variable mThreadDone_empty_cond;
  condition_variable mThreadDone_full_cond;
  condition_variable mThreadDone_initialized_cond;

  // status of the waiting list
  int mThreadDone_size; // size of the list
  int mThreadDone_pos_front; // we run jobs at the front
  int mThreadDone_pos_back; // we add new jobs at the back
  mThreadCircular<T> mThreadDone_loop;
  int mThreadDone_last_done;

  bool mThreadDone_init;
  bool mThreadDone_isrun;
};

template <typename T>
void mThreadDone<T>::stop()
{
  try {
    unique_lock<mutex> lk(mThreadDone_onebyone);
    mThreadDone_isrun = false;
    mThreadDone_empty_cond.notify_all();
  } catch (const std::exception& e) {
    std::cout << "mThreadDone<T>::stop()" << e.what();
  }
}

// initialization of the waiting list
template<typename T>
mThreadDone<T>::mThreadDone()
{
  try
  {
    mThreadDone_isrun = false;
    mThreadDone_init = false;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : mThreadDone<T>::mThreadDone()" << endl;
  }
}

// initialization of the waiting list
template<typename T>
mThreadDone<T>::mThreadDone(int size) : mThreadDone_thread( &mThreadDone<T>::run, this)
{
  try
  {
    unique_lock<mutex> lk(mThreadDone_onebyone);
    mThreadDone_isrun = true;
    mThreadDone_size = size;
    mThreadDone_pos_front = -1;
    mThreadDone_pos_back = -1;
    mThreadDone_loop.init(size);
    mThreadDone_iteration = 0;
    mThreadDone_iteration_done = 0;
    mThreadDone_last_done = -1;
    mThreadDone_init = true;
    mThreadDone_initialized_cond.notify_all();
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : mThreadDone<T>::mThreadDone(int size)" << endl;
  }
}

template<typename T>
mThreadDone<T>::~mThreadDone()
{
  try
  {
    unique_lock<mutex> full(mThreadDone_full);
    if(mThreadDone_init)
    {
      while(true_size() >= 1)
        mThreadDone_empty_cond.wait(full);
      if(mThreadDone_thread.joinable())
        mThreadDone_thread.join();
    }
    mThreadDone_init = false;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : mThreadDone<T>::~mThreadDone()" << endl;
  }
}

template <typename T>
void mThreadDone<T>::run()
{
  try
  {
    unique_lock<mutex> initialized(mThreadDone_initialized);
    while(!mThreadDone_init)
        mThreadDone_initialized_cond.wait(initialized);

    T* mThreadDone_task = nullptr;
    bool running = true;
    do
    {
      mThreadDone_task = nullptr;
      get(&mThreadDone_task); // Thread_todo->get() is supposed to block until todo is not empty
      if (mThreadDone_task != nullptr) // if todo is empty mThreadDone_task == nullptr
      {
        mThreadDone_task->done();
        delete mThreadDone_task;
      }
      else
        running = false;
    }
    while(running);
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : void mThreadDone<T>::thread_run()" << endl;
  }
}

// we add jobs at the waiting list
template<typename T>
void mThreadDone<T>::add(T* x, int n)
{
  try
  {
    if(mThreadDone_init)
    {
      unique_lock<mutex> full(mThreadDone_full);
      while(!can_add(n))
        mThreadDone_full_cond.wait(full); // if the list is full we wait
      push_back(x, n); // when there is room we add the job
    }
    else
      throw logic_error("Done list not initialized");
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : void mThreadDone<T>::add(T const & x)" << endl;
    exit(-1);
  }
}

// we load the next jobs
template<typename T>
inline void mThreadDone<T>::get(T** x)
{
  try
  {
    if(mThreadDone_init)
    {
      unique_lock<mutex> empty(mThreadDone_empty);
      while(!can_get())
        mThreadDone_empty_cond.wait(empty); // if the list in empty we wait
      *x = pop_front(); // when there is a job we get it
    }
    else
      throw logic_error("Done list not initialized");
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : T mThreadDone<T>::get()" << endl;
    exit(-1);
  }
}

// adding a job to the list
template<typename T>
inline void mThreadDone<T>::push_back(T* x, int number)
{
  try
  {
    unique_lock<mutex> lk(mThreadDone_onebyone); // we add one job at the time
    mThreadDone_loop.add(x, number);
    mThreadDone_empty_cond.notify_one(); // we signal one thread waiting for a job that the list is not empty anymore (in the get() function)
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : T mThreadDone<T>::push_back()" << endl;
    exit(-1);
  }
}

// getting a job from the list
template<typename T>
inline T* mThreadDone<T>::pop_front()
{
  try
  {
    unique_lock<mutex> lk(mThreadDone_onebyone); // we get one job at the time
    T* value = nullptr;
    value = mThreadDone_loop.pop();

    if(!mThreadDone_isrun) // if we get the stop signal there will be no new job added so their is no reason for the thread to wait for the list to fill up
    {
      mThreadDone_empty_cond.notify_all();
    }
    else // else we allow only one thread to continue
      mThreadDone_full_cond.notify_one();
    return value;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : T mThreadDone<T>::pop_front()" << endl;
    exit(-1);
  }
}

template<typename T>
inline int mThreadDone<T>::can_add(int number)
{
  try
  {
    unique_lock<mutex> lk(mThreadDone_onebyone);
    if(mThreadDone_isrun) // if we didn't get the stop signal we proceed
      return mThreadDone_loop.can_add(number);
    else // if we got the stop signal we can finish every thing even if the list is full
      throw logic_error("try to add task after the stop signal");
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : int mThreadDone<T>::can_add(int number)" << endl;
    exit(-1);
  }
}

template<typename T>
inline int mThreadDone<T>::can_get()
{
  try
  {
    unique_lock<mutex> lk(mThreadDone_onebyone);
    if(mThreadDone_isrun)
      return mThreadDone_loop.can_get();
    return true;
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : int mThreadDone<T>::can_get()" << endl;
    exit(-1);
  }
}

template<typename T>
inline int mThreadDone<T>::true_size()
{
  try
  {
    unique_lock<mutex> lk(mThreadDone_onebyone);
    return mThreadDone_loop.size();
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : int mThreadDone<T>::true_size()" << endl;
    exit(-1);
  }
}

template <typename T>
void mThreadDone<T>::join()
{
  try
  {
    unique_lock<mutex> full(mThreadDone_full);
    if(mThreadDone_init)
    {
      while(true_size() >= 1)
        mThreadDone_empty_cond.wait(full);
      if(mThreadDone_thread.joinable())
        mThreadDone_thread.join();
    }
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : int mThreadDone<T>::join()" << endl;
    exit(-1);
  }
}

template<typename T>
string mThreadDone<T>::output()
{
  try
  {
    return "D : "+to_string(mThreadDone_loop.size());
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : int mThreadDone<T>::output()" << endl;
    exit(-1);
  }
}

#endif

