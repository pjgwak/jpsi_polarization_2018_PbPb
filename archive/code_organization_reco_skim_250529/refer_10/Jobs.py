from time import sleep
import subprocess

class Jobs(object):
  def __init__(self, eventLoops):
    """ Constructor. Generates paralle jobs for the provided eventloops
    """
    self.jobNames = []
    for eventLoop in eventLoops:
      self.createExecutable(eventLoop)

  def createExecutable(self, eventLoop):
    """ for a given loop, generate stand-alone python script 
    """
    jobName = 'job.{}.py'.format(eventLoop.name)

    # open py file for writing
    with open(jobName, 'w') as f:
      # creat common steering script preamble
      f.write('from ROOT import gSystem, TH1\n')
      f.write('TH1.AddDirectory(False)\n')
      f.write('gSystem.Load("Analysis")\n')
      f.write('from Samples import *\n')
      f.write('from Algorithms import *\n')
      
      # create EventLoop class
      className = eventLoop.__class__.__name__
      f.write('eventLoop = {}()\n'.format(className))
      
      # create algorithm classes
      f.write('algs = []\n')
      for alg in eventLoop.algs:
        algClassName = alg.__class__.__name__  
        f.write('algs += [ {}() ]\n'.format(algClassName))
        f.write('eventLoop.addAlgorithms( algs) \n')
      
      # execute and save calls
      f.write('eventLoop.execute()\n')
      f.write('eventLoop.save()\n')
      f.write('print("all ok")\n')
      
      # save the jobs file name for further use
      print('Generated job executtable {}'.format(jobName))
      self.jobNames += [ jobName ]

  def execute(self, nProcesses=5, niceness=10):
      """ Executes paralle jobs
          nProcess - mas number of parallel processs
          niceness - how nice we are to other users of the computer
      """
      jobsToSubmit = self.jobNames[:] # copy the list of jobs
      runningProcess = []
      msg = ''
      while True:
        # printout
        newMsg = '{} jobs out of {} running...'.format(len(runningProcess), len(self.jobNames))
        
        if newMsg!=msg:
          print(newMsg)
          msg = newMsg
        
        # execute new job. Will only happen if number of running jobs is less than nProcess
        if len(jobsToSubmit) > 0 and len(runningProcess) < nProcesses:
          command = 'python3 {} > {}.log 2>&1'.format(jobsToSubmit[0], jobsToSubmit[0])
          print('Submitting job "{}"'.format(jobsToSubmit[0]))
          proc = subprocess.Popen('nice -n {} {}'.format(niceness,command), shell=True)
          runningProcess += [proc]
          jobsToSubmit.remove(jobsToSubmit[0])

        # check how many jobs running
        for proc in runningProcess[:]:
          if proc.poll() != None:
            # this process has terminated. Remove the running
            runningProcess.remove(proc)
        
        # exit the loop if all is submitted/done
        if len(jobsToSubmit) + len(runningProcess)==0:
          print('Done')
          break

        # wait if not in the submission stage
        if len(jobsToSubmit)==0 or len(runningProcess) >= nProcesses:
          sleep(1)