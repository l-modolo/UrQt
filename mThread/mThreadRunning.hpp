/*
htdetect horizontal transfert detection tools
Copyright (C) 2011  Laurent Modolo

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

#ifndef DEF_mThreadRunning
#define DEF_mThreadRunning

#include <future>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include "mThreadWaiting.hpp"
#include "mThreadDone.hpp"
#include <exception>
#include <stdexcept>
#include <iostream>

using namespace std;

template <typename T>
class mThreadRunning
{
  public:
  mThreadRunning(mThreadWaiting<T>* todo, mThreadDone<T>* done);
  ~mThreadRunning();

  bool get_run();
  void set_run();
  void join();
  bool joinable();

  protected:
  void run(mThreadWaiting<T>* todo, mThreadDone<T>* done);
  T* mThread_task;
  int mThread_task_number;

  mThreadWaiting<T>* mThread_todo;
  mThreadDone<T>* mThread_done;
  thread mThread_thread;
};

template <typename T>
mThreadRunning<T>::mThreadRunning() : mThread_thread(&mThreadRunning<T>::run, this, mThreadWaiting<T>* todo, mThreadDone<T>* done)
{
  // if (todo != nullptr)
  // {
  // 	mThread_todo = todo;
  // 	mThread_done = done;
  // 	mThread_task = nullptr;
  // 	mThread_task_number = -1;
  // 	// try
  // 	// {
  // 	// 	mThread_thread = thread( &mThreadRunning<T>::run, this);
  // 	// 	cout << "New thread" << endl;
  // 	// }
  // 	// catch(exception const& e)
  // 	// {
  // 	// 	cerr << "ERROR : " << e.what() << " mThreadRunning<T>::mThreadRunning(mThreadWaiting<T>* todo, mThreadDone<T>* done)" << endl;
  // 	// }
  // }
  // else
  // 	throw logic_error("mThreadRunning need a mThreadWaiting and a mThreadDone list");
}

template <typename T>
mThreadRunning<T>::~mThreadRunning()
{
  try
  {
    mThread_thread.join();
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : void mThreadRunning<T>::~mThreadRunning()" << endl;
  }
}

template <typename T>
void mThreadRunning<T>::run(mThreadWaiting<T>* todo, mThreadDone<T>* done)
{
  try
  {
    mThread_todo = todo;
    mThread_done = done;
    mThread_task = nullptr;
    mThread_task_number = -1;
    do
    {
      mThread_task = nullptr;
      cout << "get" << endl;
      mThread_todo->get(&mThread_task, &mThread_task_number); // Thread_todo->get() is supposed to block until todo is not empty
      if (mThread_task != nullptr) // if todo is empty mThread_task == nullptr
      {
        mThread_task->run();
        if(mThread_done != nullptr)
          mThread_done->add(mThread_task, mThread_task_number); // we add the object mThread_task to the done list to execute the 'mThread_task.done()' and its position in the list
                                      // Thread_done->add() is supposed to block until done is not full
      }
    }
    while(mThread_task != nullptr); // if todo is empty x == nullptr
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : void mThreadRunning<T>::run()" << endl;
  }
}

template <typename T>
void mThreadRunning<T>::join()
{
  try
  {
    mThread_thread.join();
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : void mThreadRunning<T>::join()" << endl;
  }
}

template <typename T>
bool mThreadRunning<T>::joinable()
{
  try
  {
    mThread_thread.joinable();
  }
  catch(exception const& e)
  {
    cerr << "ERROR : " << e.what() << " in : void mThreadRunning<T>::join()" << endl;
  }
}

#endif
