#!/usr/bin/python
import multiprocessing

def worker(num,l):
    print "worker",l(num)
    return


if __name__ == "__main__":
    jobs = []
    i=0
    l = lambda x:x**2
    for i in range(5):
        p = multiprocessing.Process(target = worker, args=(i,l))
        p.start()
        i+=1
