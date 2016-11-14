#!/usr/bin/python2

from os import system

def kure_job(queue, mem, job_name, job):
    job ='bsub -q {} -M {} -J {} {}'.format(
        queue,
        mem,
        job_name,
        job)
    print('Submitted Job = {}'.format(job))
    system('{}'.format(job))


def killdevil_job(queue, mem, job_name, job):
    system('bsub -q {} -M {} -J {} {}'.format(
        queue,
        mem,
        job_name,
        job))


