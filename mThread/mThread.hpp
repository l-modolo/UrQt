/*
polyNtrimmer poly nucleotide trimming tool
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

#ifndef DEF_mThread
#define DEF_mThread

#include <exception>
#include <stdexcept>
#include <iostream>
#include <future>
#include <thread>
#include <mutex>
#include <vector>
#include <condition_variable>
#include "mThreadWaiting.hpp"
#include "mThreadDone.hpp"


using namespace std;

template <typename T>
class mThread
{
	public:
	mThread(int number);
	mThread(int number, bool done_action);
	
	void stop();
	void add(T*);
	
	protected:
	void run(mThreadWaiting<T>* todo, mThreadDone<T>* done);

	bool mThread_done_task;
	mThreadWaiting<T> mThread_waiting;
	mThreadDone<T> mThread_done;
	vector< thread > mThread_running;
};

// we create a list of number thread ready to run and number*10 thread_waiting to
// pileup the jobs to do by those threads
template <typename T>
mThread<T>::mThread(int number) : mThread_waiting(number*10), mThread_done()
{
	mThread_done_task = false;
	for(int i = 0; i < number; i++)
		mThread_running.push_back( thread(&mThread<T>::run, this, &mThread_waiting, &mThread_done) );
}

// we create a list of number thread ready to run and number*10 thread_waiting to
// pileup the jobs to do by those threads
template <typename T>
mThread<T>::mThread(int number, bool done_action) : mThread_waiting(number*1000), mThread_done(number*4000)
{
	mThread_done_task = true;
	for(int i = 0; i < number; i++)
		mThread_running.push_back( thread(&mThread<T>::run, this, &mThread_waiting, &mThread_done) );
}

template <typename T>
void mThread<T>::run(mThreadWaiting<T>* todo, mThreadDone<T>* done)
{
	try
	{
		T* mThread_task = nullptr;
		int mThread_task_number = -1;
		bool running = true;
		do
		{
			mThread_task = nullptr;
			mThread_task_number = -1;
			todo->get(&mThread_task, &mThread_task_number); // Thread_todo->get() is supposed to block until todo is not empty

			if (mThread_task != nullptr) // if todo is empty mThread_task == nullptr
			{
				mThread_task->run();
				if(mThread_done_task)
					done->add(mThread_task, mThread_task_number); // we add the object mThread_task to the done list to execute the 'mThread_task.done()' and its position in the list
				else
					delete mThread_task;
			}
			else
				running = false;
		}
		while(running);
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void mThreadRunning<T>::thread_run()" << endl;
	}
}

// add a jobs to the job waiting list
template <typename T>
void mThread<T>::add(T* x)
{
	mThread_waiting.add(x);
}

// wait until all the jobs are finish
template <typename T>
void mThread<T>::stop()
{
	// we signal the jobs waiting that we wont add new jobs
	mThread_waiting.set_stop(true);
	// we wait for each jobs to finish
	for(int i = 0; i < mThread_running.size(); i++)
			(mThread_running.at(i)).join();
	mThread_done.set_stop(true);
	mThread_done.join();
}

#endif

