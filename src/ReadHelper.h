#ifndef ReadHelper_h
#define ReadHelper_h

#include <vector>
#include <string>
#include <fstream>
#include <atomic>
#include <thread>
#include <iostream>
#include <concurrentqueue.h> //https://github.com/cameron314/concurrentqueue

template <typename F>
void iterateReadsMultithreaded(const std::vector<std::string>& files, const size_t numThreads, F readCallback)
{
	std::atomic<bool> readDone;
	readDone = false;
	std::vector<std::thread> threads;
	moodycamel::ConcurrentQueue<std::shared_ptr<FastQ>> sequenceQueue;
	for (size_t i = 0; i < numThreads; i++)
	{
		threads.emplace_back([&readDone, &sequenceQueue, readCallback, i]()
		{
			while (true)
			{
				std::shared_ptr<FastQ> read;
				if (!sequenceQueue.try_dequeue(read))
				{
					bool tryBreaking = readDone;
					if (!sequenceQueue.try_dequeue(read))
					{
						if (tryBreaking) return;
						std::this_thread::sleep_for(std::chrono::milliseconds(10));
						continue;
					}
				}
				assert(read != nullptr);
				readCallback(i, *read);
			}
		});
	}
	for (const std::string& filename : files)
	{
		std::cerr << "Reading sequences from " << filename << std::endl;
		FastQ::streamFastqFromFile(filename, false, [&sequenceQueue](FastQ& read)
		{
			std::shared_ptr<FastQ> ptr = std::make_shared<FastQ>();
			std::swap(*ptr, read);
			bool queued = sequenceQueue.try_enqueue(ptr);
			if (queued) return;
			size_t triedSleeping = 0;
			while (triedSleeping < 1000)
			{
				std::this_thread::sleep_for(std::chrono::milliseconds(10));
				queued = sequenceQueue.try_enqueue(ptr);
				if (queued) return;
				triedSleeping += 1;
			}
			sequenceQueue.enqueue(ptr);
		});
	}
	readDone = true;
	for (size_t i = 0; i < threads.size(); i++)
	{
		threads[i].join();
	}
}

#endif
