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

#ifndef DEF_mThreadWaiting
#define DEF_mThreadWaiting

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
class mThreadWaiting
{
	public:
	mThreadWaiting(int size);
	~mThreadWaiting();
	
	void add(T* x);
	void get(T** x, int* number);
	
	bool stop();
	bool stop_all();
	void set_stop(bool isrun);
	void set_isrun(bool isrun);
	bool isrun();

	string output();
	
	protected:
	
	void push_back(T* x);
	T* pop_front();
	
	int size();
	int true_size();
	int mThread_iteration;
	
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
	deque<T*> mThread_waiting; // mThread_iteration of thread waiting
	
	bool mThread_isrun;
	bool mThread_stop_signal;
};

template <typename T>
bool mThreadWaiting<T>::isrun()
{
	unique_lock<mutex> lk(mThread_onebyone);
	return mThread_isrun;
}

template <typename T>
void mThreadWaiting<T>::set_isrun(bool isrun)
{
	unique_lock<mutex> lk(mThread_onebyone);
	mThread_isrun = isrun;
	if(!isrun)
		mThread_empty_cond.notify_one();
}

template <typename T>
bool mThreadWaiting<T>::stop()
{
	unique_lock<mutex> lk(mThread_onebyone);
	return mThread_stop_signal;
}

template <typename T>
bool mThreadWaiting<T>::stop_all()
{
	mThread_full_cond.notify_all();
	mThread_empty_cond.notify_all();
}

template <typename T>
void mThreadWaiting<T>::set_stop(bool isrun)
{
	unique_lock<mutex> lk(mThread_onebyone);
	mThread_stop_signal = isrun;
	mThread_empty_cond.notify_one();
}


// initialisation of the wainting list
template<typename T>
mThreadWaiting<T>::mThreadWaiting(int size)
{
	try
	{
		mThread_isrun = true;
		mThread_stop_signal = false;
		mThread_size = size;
		mThread_pos_front = -1;
		mThread_pos_back = -1;
		mThread_waiting.clear();
		mThread_iteration = -1;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : mThreadWaiting<T>::mThreadWaiting(int size)" << endl;
	}
}

template<typename T>
mThreadWaiting<T>::~mThreadWaiting()
{
	unique_lock<mutex> full(mThread_full);
	while(true_size() >= 1)
		mThread_empty_cond.wait(full);
}

// we add jobs at the waiting list
template<typename T>
void mThreadWaiting<T>::add(T* x)
{
	try
	{
		unique_lock<mutex> full(mThread_full);
		while(size() >= mThread_size)
			mThread_full_cond.wait(full); // if the list is full we wait
		push_back(x); // when there is room we add the job
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : void mThreadWaiting<T>::add(T const & x)" << endl;
		exit(-1);
	}
}

// we load the next jobs
template<typename T>
void mThreadWaiting<T>::get(T** x, int* number)
{
	try
	{
		unique_lock<mutex> empty(mThread_empty);
		while(size() <= 1 && !mThread_stop_signal)
			mThread_empty_cond.wait(empty); // if the list in empty we wait
		*x = pop_front(); // when there is a job we get it
		*number = mThread_iteration;
	}
	catch(exception const& e)
	{
		cerr << "ERROR : " << e.what() << " in : T mThreadWaiting<T>::get()" << endl;
		exit(-1);
	}
}

// adding a job to the list
template<typename T>
void mThreadWaiting<T>::push_back(T* x)
{
	unique_lock<mutex> lk(mThread_onebyone); // we add one job at the time
	mThread_waiting.push_back(x); // add job
	mThread_empty_cond.notify_one(); // we signal one thread waiting for a job that the list is not empty anymore (in the get() function)
}

// getting a job from the list
template<typename T>
T* mThreadWaiting<T>::pop_front()
{
	try
	{
		unique_lock<mutex> lk(mThread_onebyone); // we get one job at the time
		T* value = nullptr;

		if(mThread_waiting.size() <= 0)
			value = nullptr; //if the list is empty we get nothing
		else
		{
			value = mThread_waiting.at(0); // create a new object from the one at the front T need to have a initialization / copy
			mThread_iteration++;
			mThread_waiting.at(0) = nullptr;
			mThread_waiting.pop_front(); // remove the element at the front and delete it
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
		cerr << "ERROR : " << e.what() << " in : T mThreadWaiting<T>::pop_front()" << endl;
		exit(-1);
	}
}

template<typename T>
int mThreadWaiting<T>::size()
{
	unique_lock<mutex> lk(mThread_onebyone);
	if(mThread_isrun)
		return mThread_waiting.size();
	else
		return mThread_size;
}

template<typename T>
int mThreadWaiting<T>::true_size()
{
	unique_lock<mutex> lk(mThread_onebyone);
	return mThread_waiting.size();
}

template<typename T>
string mThreadWaiting<T>::output()
{
		return "W : "+to_string(mThread_waiting.size());
}

#endif

