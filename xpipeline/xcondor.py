from glue import pipeline
import os
# -------------------------------------------------------------------------
#      Define special job classes.
# -------------------------------------------------------------------------

class XsearchJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
    """
    An x search job
    """
    def __init__(self,cp):
        """
        cp = ConfigParser object from which options are read.
        """
        # ---- Get path to executable from parameters file.
        # self.__executable = cp.get('condor','xsearch')
        # ---- Get path to executable.
        os.system('which xpipeline-analysis > path_file.txt')
        f = open('path_file.txt','r')
        xpipeline_analysisstr = f.read()
        f.close()
        os.system('rm path_file.txt')
        self.__executable = xpipeline_analysisstr
        # ---- Get condor universe from parameters file.
        self.__universe = cp.get('condor','universe')
        pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
        pipeline.AnalysisJob.__init__(self,cp,False)
        self.__param_file = None

        # ---- Add required environment variables.
        self.add_condor_cmd('getenv',"true")

        # ----Add Accounting Group Flag
        grouptag = 'ligo.' + \
                cp.get('condor','ProdDevSim') + '.' +\
                cp.get('condor','Era') + '.' +\
                cp.get('condor','Group') + '.' +\
                cp.get('condor','SearchType')
        self.add_condor_cmd('accounting_group',grouptag)

        # ----Add Username Flag
        self.add_condor_cmd('accounting_group_user', cp.get('condor','UserName'))

        # ---- Add priority specification
        self.add_condor_cmd('priority', cp.get('condor', 'condorPriority'))

        # --- add minimal memory needed
        self.add_condor_cmd('request_memory', cp.get('condor', 'minimalSearchJobCPUMem'))

        # ---- Path and file names for standard out, standard error for this job.
        self.set_stdout_file('logs/xsearch-$(cluster)-$(process).out')
        self.set_stderr_file('logs/xsearch-$(cluster)-$(process).err')

        # ---- Name of condor job submission file to be written.
        self.set_sub_file('xsearch.sub')

class XsearchGPUJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
    """
    An x search job
    """
    def __init__(self,cp):
        """
        cp = ConfigParser object from which options are read.
        """
        # FOR GPUS

        # ---- Get path to executable from parameters file.
        # self.__executable = cp.get('condor','xsearch')
        # ---- Get path to executable.
        os.system('which xpipeline-analysis > path_file.txt')
        f = open('path_file.txt','r')
        xpipeline_analysisstr = f.read()
        f.close()
        os.system('rm path_file.txt')
        self.__executable = xpipeline_analysisstr
        # ---- Get condor universe from parameters file.
        self.__universe = cp.get('condor','universe')
        pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
        pipeline.AnalysisJob.__init__(self,cp,False)
        self.__param_file = None

        # ---- Add required environment variables.
        self.add_condor_cmd('environment',"USER=$ENV(USER);HOME=$ENV(HOME);" \
            "LD_LIBRARY_PATH=$ENV(LD_LIBRARY_PATH);" \
            "XPIPE_INSTALL_BIN=$ENV(XPIPE_INSTALL_BIN);" \
            "PATH=/usr/bin:/bin")

        # ----Add Accounting Group Flag
        grouptag = 'ligo.' + \
                cp.get('condor','ProdDevSim') + '.' +\
                cp.get('condor','Era') + '.' +\
                cp.get('condor','Group') + '.' +\
                cp.get('condor','SearchType')
        self.add_condor_cmd('accounting_group',grouptag)

        # ----Add Username Flag
        self.add_condor_cmd('accounting_group_user',cp.get('condor','UserName'))

        # ---- Add priority specification
        self.add_condor_cmd('priority', cp.get('condor', 'condorPriority'))

        # --- add minimal memory needed
        self.add_condor_cmd('request_memory', cp.get('condor', 'minimalSearchJobGPUMem'))

        # ---- Path and file names for standard out, standard error for this job.
        self.set_stdout_file('logs/xsearch-$(cluster)-$(process).out')
        self.set_stderr_file('logs/xsearch-$(cluster)-$(process).err')

        # If on Atlas, use getenv=true to pass variables
        #if atlasFlag:
        #   self.add_condor_cmd('getenv',"true")

        self.add_condor_cmd('Requirements','TARGET.WantGPU =?= True')
        self.add_condor_cmd('+WantGPU','True')

        # ---- Name of condor job submission file to be written.
        self.set_sub_file('xsearch_gpu.sub')


class XsearchNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
    """
    xsearch node
    """
    def __init__(self,job):
        """
        job = A CondorDAGJob.
        """
        pipeline.CondorDAGNode.__init__(self,job)
        pipeline.AnalysisNode.__init__(self)
        self.__x_jobnum = None
        self.__x_injnum = None

    # ---- Set parameters file.
    def set_param_file(self,path):
        self.add_var_arg(path)
        self.__param_file = path

    def get_param_file(self):
        return self.__param_file

    def set_x_jobnum(self,n):
        self.add_var_arg(str(n))
        self.__x_jobnum = n

    def get_x_jobnum(self):
        return self.__x_jobnum

    def set_output_dir(self,path):
        self.add_var_arg(path)
        self.__output_dir = path

    def get_output_dir(self,path):
        return self.__output_dir

    def set_x_injnum(self,n):
        self.add_var_arg(n)
        self.__x_injnum = n

    def get_x_injnum(self):
        return self.__x_injnum

class XmergeJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
    """
    An x merge job
    """
    def __init__(self,cp):
        """
        cp = ConfigParser object from which options are read.
        """
        # ---- Get path to executable.
        os.system('which xmergegrbresults > path_file.txt')
        f = open('path_file.txt','r')
        xmergestr = f.read()
        f.close()
        os.system('rm path_file.txt')
        self.__executable = xmergestr
        # ---- Get condor universe from parameters file.
        self.__universe = cp.get('condor','universe')
        pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
        pipeline.AnalysisJob.__init__(self,cp,False)
        # ---- Add name of 'output' directory as first argument.
        self.add_arg('output')
        # ---- Add required environment variables.
        self.add_condor_cmd('environment',"USER=$ENV(USER);HOME=$ENV(HOME);" \
            "LD_LIBRARY_PATH=$ENV(LD_LIBRARY_PATH)")

        # ----Add Accounting Group Flag
        grouptag = 'ligo.' + \
                cp.get('condor','ProdDevSim') + '.' +\
                cp.get('condor','Era') + '.' +\
                cp.get('condor','Group') + '.' +\
                cp.get('condor','SearchType')
        self.add_condor_cmd('accounting_group',grouptag)

        # ----Add Username Flag
        self.add_condor_cmd('accounting_group_user',cp.get('condor','UserName'))

        # ---- Add priority specification
        self.add_condor_cmd('priority', cp.get('condor', 'condorPriority'))

        # --- add minimal memory needed
        self.add_condor_cmd('request_memory', cp.get('condor', 'minimalMergeJobMem'))

        # ---- Path and file names for standard out, standard error for this job.
        self.set_stdout_file('logs/xmerge-$(cluster)-$(process).out')
        self.set_stderr_file('logs/xmerge-$(cluster)-$(process).err')

        # If on Atlas
        #if atlasFlag:
        #   self.add_condor_cmd('getenv',"true")

        # ---- Name of condor job submission file to be written.
        self.set_sub_file('xmerge.sub')

class XmergeNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
    """
    merge results that were cut into blocks of less than maxInjNum jobs
    """
    def __init__(self,job):
        pipeline.CondorDAGNode.__init__(self,job)
        pipeline.AnalysisNode.__init__(self)

    def set_dir_prefix(self,path):
        self.add_var_arg(path)
        self.__dir_prefix = path

    def get_dir_prefix(self):
        return self.__dir_prefix

    def set_sn_flag(self,path):
        self.add_var_arg(path)
        self.__sn_flag = path

    def get_sn_flag(self):
        return self.__sn_flag


class XmergeClusteredJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  An x merge clustered job
  """
  def __init__(self,cp,nDet):
    # ---- Get path to executable.
    """
    cp = ConfigParser object from which options are read.
    """
    if nDet==3:
      os.system('which xmergegrbclusteredresults > path_file.txt')
    else:
      os.system('which xmergegrbclusteredresultsTwoDets > path_file.txt')
    f = open('path_file.txt','r')
    xmergestr = f.read()
    f.close()
    os.system('rm path_file.txt')
    self.__executable = xmergestr
    # ---- Get condor universe from parameters file.
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,False)
    # ---- Add name of 'output' directory as first argument.
    self.add_arg('output_clustered')
    # ---- Add required environment variables.
    self.add_condor_cmd('environment',"USER=$ENV(USER);HOME=$ENV(HOME);" \
        "LD_LIBRARY_PATH=$ENV(LD_LIBRARY_PATH)")

    # ----Add Accounting Group Flag
    grouptag = 'ligo.' + \
            cp.get('condor','ProdDevSim') + '.' +\
            cp.get('condor','Era') + '.' +\
            cp.get('condor','Group') + '.' +\
            cp.get('condor','SearchType')
    self.add_condor_cmd('accounting_group',grouptag)

    # ----Add Username Flag
    self.add_condor_cmd('accounting_group_user',cp.get('condor','UserName'))

    # ---- Add priority specification.
    self.add_condor_cmd('priority', cp.get('condor', 'condorPriority'))

    # ---- Add minimal memory request. Note that with current default minimalMem
    #      we never access the "else" statement below.
    self.add_condor_cmd('request_memory', cp.get('condor', 'minimalMergeClusterJobMem'))

    # ---- Path and file names for standard out, standard error for this job.
    self.set_stdout_file('logs/xmergeclustered-$(cluster)-$(process).out')
    self.set_stderr_file('logs/xmergeclustered-$(cluster)-$(process).err')

    # If on Atlas
    #if atlasFlag:
    # self.add_condor_cmd('getenv',"true")

    # ---- Name of condor job submission file to be written.
    self.set_sub_file('xmergeclustered.sub')

class XmergeClusteredNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  cluster and merge results that were cut into blocks of less than maxInjNum jobs
  """
  def __init__(self,job):
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)

  def set_dir_prefix(self,path):
    self.add_var_arg(path)
    self.__dir_prefix = path

  def get_dir_prefix(self):
    return self.__dir_prefix

  def set_sn_flag(self,path):
    self.add_var_arg(path)
    self.__sn_flag = path

  def get_sn_flag(self):
    return self.__sn_flag

  def set_sc_flag(self,path):
    self.add_var_arg(path)
    self.__sc_flag = path

  def get_sc_flag(self):
    return self.__sc_flag

class XtmvaJob(pipeline.CondorDAGJob, pipeline.AnalysisJob):
  """
  An xtmva classification job.
  """
  def __init__(self,cp,nDet):
    """
    cp = ConfigParser object from which options are read.
    """
    # ---- Get path to executable.
    self.__executable = os.getcwd() + "/xtmvapy/xtmva.py \n"
    # ---- Get condor universe from parameters file.
    self.__universe = cp.get('condor','universe')
    pipeline.CondorDAGJob.__init__(self,self.__universe,self.__executable)
    pipeline.AnalysisJob.__init__(self,cp,False)
    # ---- Add name of 'output' directory as first argument.
    cwdstr_split = os.path.split(os.getcwd())
    fdir = cwdstr_split[0]
    grbdir = cwdstr_split[1]
    argumentstr = " --fdir=" + fdir + " --grbname=" + grbdir + " --nifo=" + str(nDet) + " --cname=xtmva.ini "
    print >> sys.stdout, argumentstr
    # ---- Add name of 'output' directory as first argument.
    self.add_arg(argumentstr)
    # ---- Define environment variables that are needed by xtmva.py. We hard-code 
    #      most of them to avoid need to create a special gen-wrapper-script*.sh 
    #      script for this one function.
    envarstr = "USER=" + os.getenv("USER") + ";" \
    + "HOME=" + os.getenv("HOME") + ";" \
    + "LD_LIBRARY_PATH=" + os.getenv("LD_LIBRARY_PATH") + ";" \
    + "XPIPE_INSTALL_BIN=$ENV(XPIPE_INSTALL_BIN);" \
    + "PATH=" + os.getenv("ROOTSYS") + "/bin:" + os.getenv("XPIPE_INSTALL_BIN") + ":/usr/bin:/bin;" \
    + "LIBPATH=" + os.getenv("LIBPATH") + ";" \
    + "DYLD_LIBRARY_PATH=" + os.getenv("DYLD_LIBRARY_PATH") + ";" \
    + "PYTHONPATH=" + os.getenv("PYTHONPATH") + ";" \
    + "SHLIB_PATH=" +  os.getenv("SHLIB_PATH") + ";" \
    + "ROOTSYS=" + os.getenv("ROOTSYS")
    # ---- Add required environment variables.
    self.add_condor_cmd('environment',envarstr)

    # ---- Add priority specification.
    self.add_condor_cmd('priority', cp.get('condor', 'condorPriority'))

    # ---- Add minimal memory request. 
    self.add_condor_cmd('request_memory', cp.get('condor', 'minimalXtmvaJobMem'))

    # ----Add Accounting Group Flag
    grouptag = 'ligo.' + \
            cp.get('condor','ProdDevSim') + '.' +\
            cp.get('condor','Era') + '.' +\
            cp.get('condor','Group') + '.' +\
            cp.get('condor','SearchType')
    self.add_condor_cmd('accounting_group',grouptag)

    # ----Add Username Flag
    self.add_condor_cmd('accounting_group_user',cp.get('condor','UserName'))

    # ---- Path and file names for standard out, standard error for this job.
    self.set_stdout_file('logs/xtmva-$(cluster)-$(process).out')
    self.set_stderr_file('logs/xtmva-$(cluster)-$(process).err')

    # If on Atlas
    #if atlasFlag:
    # self.add_condor_cmd('getenv',"true")

    # ---- Name of condor job submission file to be written.
    self.set_sub_file('xtmva.sub')

class XtmvaNode(pipeline.CondorDAGNode, pipeline.AnalysisNode):
  """
  cluster and merge results that were cut into blocks of less than maxInjNum jobs
  """
  def __init__(self,job):
    pipeline.CondorDAGNode.__init__(self,job)
    pipeline.AnalysisNode.__init__(self)
