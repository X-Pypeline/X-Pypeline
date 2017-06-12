def parse(outputType, cp):
    # ---- xdetection is not yet set up to read and use this clustering parameters
    #      file.  Therefore, this section of setUpJobs is commented out until
    #      xdetection is able to use it.

    # # ---- Note: If outputType is clusters, then the parameters ini file must
    # #      contain a section [extraction] which specifies the manner in which
    # #      clustering is to be performed.
    # outputType=cp.get('parameters','outputType')
    # if outputType == 'clusters' :
    #
    #     # ---- Status message.
    #     print >> sys.stdout, "Writing Matlab-formatted clustering file ..."
    #     extractParamFileName = cp.get('parameters','extractParamFileName')
    #
    #     # ---- Parameters file for on-source analysis.
    #     f=open(extractParamFileName, 'w')
    #
    #     # ---- Write all options available in the [extraction] section to a file.
    #     extractParams = cp.options('extraction')
    #     for i in range(0,len(extractParams)) :
    #         value = cp.get('extraction',extractParams[i])
    #         f.write(extractParams[i] + ':' + value + '\n')
    #     f.close()
    #
    # # ---- Status message.
    # print >> sys.stdout, "... finished writing clustering file.     "
    # print >> sys.stdout

    if outputType == 'seedless':
	 # ---- Status message.
	 print >> sys.stdout, "Writing Matlab-formatted clustering file ..."
	 doGPUOnSource = cp.get('seedless','doGPUOnSource')
	 doGPUOffSource = cp.get('seedless','doGPUOffSource')
	 doGPUSimulation = cp.get('seedless','doGPUSimulation')
	 doGPUUL = cp.get('seedless','doGPUUL')
	 doParallelOnSource = cp.get('seedless','doParallelOnSource')
	 doParallelOffSource = cp.get('seedless','doParallelOffSource')
	 doParallelSimulation = cp.get('seedless','doParallelSimulation')
	 doParallelUL = cp.get('seedless','doParallelUL')
	 T = cp.get('seedless','T')
	 F = cp.get('seedless','F')
	 mindur = cp.get('seedless','mindur')
	 doPCA = cp.get('seedless','doPCA')
	 pca_catalogfile = cp.get('seedless','pca_catalogfile')
	 pca_type = cp.get('seedless','pca_type')
	 doCBC = cp.get('seedless','doCBC')
	 doECBC = cp.get('seedless','doECBC')
	 doBezier = cp.get('seedless','doBezier')
	 doExponential = cp.get('seedless','doExponential')
	 doRModes = cp.get('seedless','doRModes')
	 norm = cp.get('seedless','norm')
	 savePlots = cp.get('seedless','savePlots')


	 f = open('input/seedless_onsource.txt','w')
	 f.write('doGPU:%s\n'%(doGPUOnSource))
	 f.write('doParallel:%s\n'%(doParallelOnSource))
	 f.write('T:%s\n'%(T))
	 f.write('F:%s\n'%(F))
	 f.write('mindur:%s\n'%(mindur))
	 f.write('doCBC:%s\n'%(doCBC))
	 f.write('doECBC:%s\n'%(doECBC))
	 f.write('doBezier:%s\n'%(doBezier))
	 f.write('doExponential:%s\n'%(doExponential))
	 f.write('doRModes:%s\n'%(doRModes))
	 f.write('doPCA:%s\n'%(doPCA))
	 f.write('pca_catalogfile:%s\n'%(pca_catalogfile))
	 f.write('pca_type:%s\n'%(pca_type))
	 f.write('norm:%s\n'%(norm))
	 f.write('savePlots:%s\n'%(savePlots))
	 f.close()

	 f = open('input/seedless_offsource.txt','w')
	 f.write('doGPU:%s\n'%(doGPUOffSource))
	 f.write('doParallel:%s\n'%(doParallelOffSource))
	 f.write('T:%s\n'%(T))
	 f.write('F:%s\n'%(F))
	 f.write('mindur:%s\n'%(mindur))
	 f.write('doCBC:%s\n'%(doCBC))
	 f.write('doECBC:%s\n'%(doECBC))
	 f.write('doBezier:%s\n'%(doBezier))
	 f.write('doExponential:%s\n'%(doExponential))
	 f.write('doRModes:%s\n'%(doRModes))
	 f.write('doPCA:%s\n'%(doPCA))
	 f.write('pca_catalogfile:%s\n'%(pca_catalogfile))
	 f.write('pca_type:%s\n'%(pca_type))
	 f.write('norm:%s\n'%(norm))
	 f.write('savePlots:%s\n'%(savePlots))
	 f.close()

	 f = open('input/seedless_simulation.txt','w')
	 f.write('doGPU:%s\n'%(doGPUSimulation))
	 f.write('doParallel:%s\n'%(doParallelSimulation))
	 f.write('T:%s\n'%(T))
	 f.write('F:%s\n'%(F))
	 f.write('mindur:%s\n'%(mindur))
	 f.write('doCBC:%s\n'%(doCBC))
	 f.write('doECBC:%s\n'%(doECBC))
	 f.write('doBezier:%s\n'%(doBezier))
	 f.write('doExponential:%s\n'%(doExponential))
	 f.write('doRModes:%s\n'%(doRModes))
	 f.write('doPCA:%s\n'%(doPCA))
	 f.write('pca_catalogfile:%s\n'%(pca_catalogfile))
	 f.write('pca_type:%s\n'%(pca_type))
	 f.write('norm:%s\n'%(norm))
	 f.write('savePlots:%s\n'%(savePlots))
	 f.close()

	 f = open('input/seedless_ul.txt','w')
	 f.write('doGPU:%s\n'%(doGPUUL))
	 f.write('doParallel:%s\n'%(doParallelUL))
	 f.write('T:%s\n'%(T))
	 f.write('F:%s\n'%(F))
	 f.write('mindur:%s\n'%(mindur))
	 f.write('doCBC:%s\n'%(doCBC))
	 f.write('doECBC:%s\n'%(doECBC))
	 f.write('doBezier:%s\n'%(doBezier))
	 f.write('doExponential:%s\n'%(doExponential))
	 f.write('doRModes:%s\n'%(doRModes))
	 f.write('doPCA:%s\n'%(doPCA))
	 f.write('pca_catalogfile:%s\n'%(pca_catalogfile))
	 f.write('pca_type:%s\n'%(pca_type))
	 f.write('norm:%s\n'%(norm))
	 f.write('savePlots:%s\n'%(savePlots))
	 f.close()
