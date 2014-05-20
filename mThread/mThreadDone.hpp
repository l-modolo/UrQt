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

#ifndef DEF_mThreadDone
#define DEF_mThreadDone

#include <future>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <deque>
#include <exception>
#include <stdexcept>
#include <iostream>
#include <string>

using namespace std;

template <typename T>
class mThreadDone
{
	public:
	mThreadDone();
	mThreadDone(int size);
	~mThreadDone();
	
	void add(T* x, int number);
	void get(T** x);
	
	bool stop();
	bool stop_all();
	void set_stop(bool isrun);
	void set_isrun(bool isrun);
	bool isrun();
	
	void join();
	
	string output();

	protected:
	void push_back(T* x, int number);
	T* pop_front();
	
	int can_add(int number);
	int can_get();
	int true_size();
	int mThread_iteration;
	int mThread_iteration_done;
	
	thread mThread_thread;
	void run();
	
	// communication between the thread
	mutex mThread_onebyone;
	mutex mThread_empty; // if the list is empty we wait new jobs
	mutex mThread_full; // if the list is full we wait before new jobs
	condition_variable mThread_empty_cond;
	condition_variable mThread_full_cond;
	
	// status of the waiting list
	int mThread_size; // size of the list
	int mThread_pos_front; // we isrun jobs at the front
	int mThread_pos_back; // we add new jobs at the back
	deque<T*> mThread_waiting; // thread waiting
	deque<int> mThread_waiting_number; // id of thread waiting
	int mThread_last_done;
	
	bool mThread_init;
	bool mThread_isrun;
	bool mThread_stop_signal;
};

template <typename T>
bool mThreadDone<T>::isrun()
{
	unique_lock<mutex> lk(mThread_onebyone);
	return mThread_isrun;
}

template <typename T>
void mThreadDone<T>::set_isrun(bool isrun)
{
	unique_lock<mutex> lk(mThread_onebyone);
	mThread_isrun = isrun;
	if(!isrun)
		mThread_empty_cond.notify_one();
}

template <typename T>
bool mThreadDone<T>::stop()
{
	unique_lock<mutex> lk(mThread_onebyone);
	return mThread_stop_signal;
}

template <typename T>
bool mThreadDone<T>::stop_all()
{
	mThread_full_cond.notify_all();
	mThread_empty_cond.notify_all();
}

template <typename T>
void mThreadDone<T>::set_stop(bool isrun)
{
	unique_lock<mutex> lk(mThread_onebyone);
	mThread_stop_signal = isrun;
	mThread_empty_cond.notify_one();
}

// initialization of the waiting list
template<typename T>
mThreadDone<T>::mThreadDone()
{
	try
	{
		unique_lock<mutex> lk(mThread_onebyone);
		mThread_isrun = false;
		mThread_stop_signal = true;
		mThread_init = false;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : mThreadDone<T>::mThreadDone()" << endl;
	}
}

// initialization of the waiting list
template<typename T>
mThreadDone<T>::mThreadDone(int size) : mThread_thread( &mThreadDone<T>::run, this)
{
	try
	{
		unique_lock<mutex> lk(mThread_onebyone);
		mThread_isrun = true;
		mThread_stop_signal = false;
		mThread_size = size;
		mThread_pos_front = -1;
		mThread_pos_back = -1;
		mThread_waiting.clear();
		mThread_iteration = 0;
		mThread_iteration_done = 0;
		mThread_last_done = -1;
		mThread_init = true;
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
		unique_lock<mutex> full(mThread_full);
		if(mThread_init)
		{
				while(true_size() >= 1)
					mThread_empty_cond.wait(full);
				// cout << "done finish" << endl;
				join();
		}
		mThread_init = false;
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
		T* mThread_task = nullptr;
		int mThread_task_number = -1;
		do
		{
			
			mThread_task = nullptr;
			get(&mThread_task); // Thread_todo->get() is supposed to block until todo is not empty
			if (mThread_task != nullptr) // if todo is empty mThread_task == nullptr
			{
				mThread_task->done();
				delete mThread_task;
			}
		}
		while(mThread_task != nullptr); // if todo is empty x == nullptr
		// cout << "writing done" << endl;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void mThreadRunning<T>::thread_run()" << endl;
	}
}

// we add jobs at the waiting list
template<typename T>
void mThreadDone<T>::add(T* x, int n)
{
	try
	{
		if(mThread_init)
		{
			int number = n;
			unique_lock<mutex> full(mThread_full);
			while(!can_add(number))
				mThread_full_cond.wait(full); // if the list is full we wait
			push_back(x, number); // when there is room we add the job
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
void mThreadDone<T>::get(T** x)
{
	try
	{
		if(mThread_init)
		{
			unique_lock<mutex> empty(mThread_empty);
			while(!can_get() && !mThread_stop_signal)
				mThread_empty_cond.wait(empty); // if the list in empty we wait
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
void mThreadDone<T>::push_back(T* x, int number)
{
	unique_lock<mutex> lk(mThread_onebyone); // we add one job at the time
	if(mThread_waiting_number.size() == 0 || number > mThread_waiting_number.back()) // if the job is the last to be executed we add it at the back
	{
		mThread_waiting.push_back(x); // add job at the back of the list
		mThread_waiting_number.push_back(number);
	}
	else // else we add it at the right position in the list
	{
		typename deque< T* >::iterator it = mThread_waiting.begin();
		deque< int >::iterator it_number = mThread_waiting_number.begin();
		while (it_number != mThread_waiting_number.end() and *it_number < number)
		{
			it++;
			it_number++;
		}
		mThread_waiting.insert(it,x); // add job just before the next job
		mThread_waiting_number.insert(it_number,number);

	}
	mThread_empty_cond.notify_one(); // we signal one thread waiting for a job that the list is not empty anymore (in the get() function)
}

// getting a job from the list
template<typename T>
T* mThreadDone<T>::pop_front()
{
	try
	{
		unique_lock<mutex> lk(mThread_onebyone); // we get one job at the time
		T* value = nullptr;

		if(mThread_waiting.size() <= 0)
			value = nullptr; //if the list is empty we get nothing
		else
		{
			value = mThread_waiting.at(0);
			mThread_waiting.at(0) = nullptr; // create a new object from the one at the front T need to have a initiator / copy
			mThread_waiting.pop_front(); // remove the element at the front and delete it
			mThread_last_done = mThread_waiting_number.at(0); // we record the number of the last job done
			mThread_waiting_number.pop_front();
		}

		if(mThread_stop_signal) // if we get the stop signal there will be no new job added so their is no reason for the thread to wait for the list to fill up
		{
			mThread_isrun = false;
			mThread_empty_cond.notify_all();
		}
		else // else we allow only one thread to continue
			mThread_full_cond.notify_one();

		return value;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : T mThreadDone<T>::pop_front()" << endl;
		exit(-1);
	}
}

template<typename T>
int mThreadDone<T>::can_add(int number)
{
	unique_lock<mutex> lk(mThread_onebyone);
	if(mThread_isrun) // if we didn't get the stop signal we proceed
	{
		if(mThread_waiting.size() < mThread_size) // if the list is not full
			return true;
		else // if the list is full we check if the number jobs number is after the max job number of the list
		{
			if(mThread_waiting_number.back() == -1)
				return true;
			if(number < mThread_waiting_number.back())
				return true;
			else
				return false;
		}
	}
	else // if we got the stop signal we can finish every thing even if the list is full
		return true;
}

template<typename T>
int mThreadDone<T>::can_get()
{
	try
	{
		unique_lock<mutex> lk(mThread_onebyone);
		if(mThread_waiting.size() > 0) // if the list is not empty
		{
			if(mThread_last_done == -1 && mThread_waiting_number.front() == 0) // if it's the first job we run it
				return true;
			else // if we already run the first job
			{
				if(mThread_last_done == mThread_waiting_number.front() - 1) // if it's the next jobs we run it
					return true;
				else // else we continue to wait for the next job
					return false;
			}
		}
		else // if the list is empty we wait
		{
			if(mThread_isrun) // if we didn't get the stop signal we proceed
				return false;
			else
				return true;
		}
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : int mThreadDone<T>::can_get()" << endl;
		exit(-1);
	}
}

template<typename T>
int mThreadDone<T>::true_size()
{
	unique_lock<mutex> lk(mThread_onebyone);
	return mThread_waiting.size();
}

template <typename T>
void mThreadDone<T>::join()
{
	if(mThread_init)
		if(mThread_thread.joinable())
			mThread_thread.join();
}

template<typename T>
string mThreadDone<T>::output()
{
	return "D : "+to_string(mThread_waiting.size());
}

#endif

