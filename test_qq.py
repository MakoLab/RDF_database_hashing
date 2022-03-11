import queue


pq = queue.PriorityQueue()
pq.put(((1, 5, 4, "ala", "beta"), "pierwszy"))
pq.put((("mala", 4, 6, 3), "drugi"))
pq.put(((1, 2, 7, "nic", "kot"), "trzeci"))

print(pq)
