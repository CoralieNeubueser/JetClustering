#!/usr/bin/env python
import os, sys, subprocess, time, commands
from optparse import OptionParser

#__________________________________________________________
def getCommandOutput(command):
    p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout,stderr = p.communicate()
    return {"stdout":stdout, "stderr":stderr, "returncode":p.returncode}

#__________________________________________________________
def SubmitToCondor(cmd,nbtrials):
    submissionStatus=0
    cmd=cmd.replace('//','/')
    for i in xrange(nbtrials):
        outputCMD = getCommandOutput(cmd)
        stderr=outputCMD["stderr"].split('\n')
        stdout=outputCMD["stdout"].split('\n')

        if len(stderr)==1 and stderr[0]=='' :
            print "------------GOOD SUB"
            submissionStatus=1
        else:
            print "++++++++++++ERROR submitting, will retry"
            print "Trial : "+str(i)+" / "+str(nbtrials)
            print "stderr : ",stderr
            print "stderr : ",len(stderr)

            time.sleep(10)


        if submissionStatus==1:
            return 1,0

        if i==nbtrials-1:
            print "failed sumbmitting after: "+str(nbtrials)+" trials, will exit"
            return 0,0

def main():
    parser = OptionParser()

    parser.add_option ('-i','--input', help='input directory containing hgcal ntuple trees (most likely on eos)',
                       dest='input',
                       default='')

    parser.add_option ('-p','--pt', help='pt',
                       dest='pt',
                       default='10')

    parser.add_option ('-q','--queue', help='lsf queue',
                       dest='queue',
                       default='1nd')

    parser.add_option ('-n', '--nev',  help='max number of event per file. default is 1000',
                       dest='nev',
                       default='1000')

    parser.add_option ('--njobs',  help='max number of jobs. Max is number of ntuple files',
                       dest='njobs',
                       default='10')

    parser.add_option ('--deltaR',  help='cone size of anti-kt',
                       dest='deltaR',
                       default='0.4')

    parser.add_option ('--etaCut',  help='cone size of anti-kt',
                       dest='etaCut',
                       default='1.3')

    parser.add_option("--algorithm", help="choose algorithm : antiKt, antiKt_cluster, or simpleCone",
                      dest="algorithm",
                      default='antiKt')

    parser.add_option("-c","--collect", help="collects jobs for given algorithm",
		      dest="collect", 
		      default='')

    parser.add_option("-a","--allPts", help="collects jobs for all pts",
		      dest="allPts", action="store_true", 
		      default=False)

    parser.add_option("-t","--topoCluster", help="topoCluster for jet building",
                      dest="topoCluster", action="store_true",
                      default=False)

    parser.add_option("--pileupNoise", help="added pileup noise before jet building",
                      dest="pileupNoise", type=int,
                      default=0)

    parser.add_option("--coneCheck", help="check simple cone sum around jet axis",
                      dest="coneCheck",
                      default="0")

    parser.add_option("--bFieldOff", help="b field is switched off",
                      dest="bFieldOff", action="store_true",
                      default=False)

    parser.add_option("--version", help="version of FCCSimJobs to find/define input/output",
                      dest="version",
                      default="v03")

    parser.add_option("--etaMax", help="etaMax",
                      dest="etaMax",
                      default=1.5)

    parser.add_option("--process", help="process",
                      dest="proc",
                      default="ljets")

    (options, args) = parser.parse_args()
    input_dir       = options.input
    algo            = options.algorithm
    deltaR = options.deltaR
    etaCut = options.etaCut
    bfield = True
    bfieldText = 'bFieldOn'

    if options.bFieldOff:
        bfield = False
        bfieldText = 'bFieldOff'
    check = "0"

    if options.coneCheck == "1":
        algo = algo + '_coneCheck'
        check = "1"

    if options.topoCluster:
        algo = algo + '_cluster'

    if options.pileupNoise:
        algo = algo + '_pileupNoise_mu' + str(options.pileupNoise)

    output_dir    = '/eos/experiment/fcc/hh/simulation/samples/'+options.version+'/physics/'+str(options.proc)+'/'+bfieldText+'/etaTo'+str(options.etaMax)+'/'+options.pt+'GeV/ana/JetClustering/'+algo+'/DeltaR_'+deltaR+'/maxEta'+etaCut+'/'
    max_events    = int(options.nev)
    queue         = options.queue
    collect       = options.collect
    condorPath    = '/afs/cern.ch/user/c/cneubuse/FCC/JetClustering/batch'

    if collect and not options.allPts:
        en = int(options.pt)
        print 'Collecting jobs for process: '+str(collect)
        input_dir = output_dir
        # /eos/experiment/fcc/hh/simulation/samples/v03/physics/'+str(options.proc)+'/bFieldOn/etaTo1.5/100GeV/ana/JetClustering/anti-kt_coneCheck/out/
        print 'Find files in '+input_dir+'/out/'
        hadd_dir = input_dir+'/out/'
        basename = str(collect)+'_'+str(options.pt)+'GeV_'+str(options.proc)
        basename = os.path.basename(basename)
        print basename
        outfile = hadd_dir + basename + '.root'
        hadd_files = hadd_dir + '*.root'
        cmd ='hadd -f -k -n 0 '+ outfile + ' ' + hadd_files
        os.system(cmd)
        sys.exit('Collection of jobs done.')
       
    if collect and options.allPts:
        energies = [20,50,100,200,500,1000,2000,5000,10000]
        outfile = '~/FCC/JetClustering/anti-kt_inBfield/'+options.version+'/DeltaR_'+deltaR+'/maxEta'+etaCut+'/all_'+collect+'.root'
        if options.bFieldOff:
            outfile = '~/FCC/JetClustering/anti-kt_noBfield/'+options.version+'/DeltaR_'+deltaR+'/maxEta'+etaCut+'/all_'+collect+'.root'
        elif options.version == 'v02_pre':
            outfile = '~/FCC/JetClustering/anti-kt_noBfield/'+options.version+'/DeltaR_'+deltaR+'/maxEta'+etaCut+'/all_'+collect+'.root'
        print outfile
        for energy in energies:
            dir_out = '/eos/experiment/fcc/hh/simulation/samples/'+options.version+'/physics/'+str(options.proc)+'/'+bfieldText+'/etaTo'+str(options.etaMax)+'/'+str(energy)+'GeV/ana/JetClustering/'+collect+'/DeltaR_'+deltaR+'/maxEta'+etaCut+'/out/'
#            if float(deltaR) == 0.4 :
#                input_dir = '/eos/experiment/fcc/hh/simulation/samples/'+options.version+'/physics/'+str(options.proc)+'/'+bfieldText+'/etaTo1.5/'+str(energy)+'GeV/ana/JetClustering/'+collect+'/out/'
#                dir_out = input_dir
            print 'directory in which to find the merged files of energy {0}GeV : {1}'.format(energy,dir_out)
            if options.bFieldOff:
                os.system('cp '+dir_out+'/'+str(collect)+'_'+str(energy)+'GeV.root ~/FCC/JetClustering/anti-kt_noBfield/'+options.version+'/DeltaR_'+deltaR+'/maxEta'+etaCut+'/'+str(collect)+'_'+str(energy)+'GeV.root')
            elif options.version == 'v02_pre':
                os.system('cp '+dir_out+'/'+str(collect)+'_'+str(energy)+'GeV.root ~/FCC/JetClustering/anti-kt_noBfield/'+options.version+'/DeltaR_'+deltaR+'/maxEta'+etaCut+'/')
            else:
                os.system('cp '+dir_out+'/'+str(collect)+'_'+str(energy)+'GeV.root ~/FCC/JetClustering/anti-kt_inBfield/'+options.version+'/DeltaR_'+deltaR+'maxEta'+etaCut+'/')
        
        if options.bFieldOff:
            hadd_files = '~/FCC/JetClustering/anti-kt_noBfield/'+options.version+'/DeltaR_'+deltaR+'/maxEta'+etaCut+'/'+str(collect)+'*GeV.root'
        elif options.version == 'v02_pre':
            hadd_files = '~/FCC/JetClustering/anti-kt_noBfield/'+options.version+'/DeltaR_'+deltaR+'/maxEta'+etaCut+'/'+str(collect)+'*GeV.root'
        else:
            hadd_files = '~/FCC/JetClustering/anti-kt_inBfield/'+options.version+'/DeltaR_'+deltaR+'/maxEta'+etaCut+'/'+str(collect)+'*GeV.root'
        os.system('hadd -f -k '+ outfile + ' ' + hadd_files)
        sys.exit('All energies have been collected, saved as: '+outfile)

    # first create output dir
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    if not os.path.exists(condorPath+'/err/'):
        os.makedirs(condorPath+'/err/')
        os.makedirs(condorPath+'/log/')
        os.makedirs(condorPath+'/out/')
    else:
       print 'Output dir: "{0}" exists.'.format(output_dir)

    cmd = "find {0} -type f -name '*.root'".format(options.input)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    lst = process.communicate()[0]
    list_of_files = lst.splitlines()
    
    # just send one job per ntup file
    nfiles = len(list_of_files)
    njobs = min(nfiles, int(options.njobs))
    
    for job in xrange(njobs):
       
       print "output dir  : ", output_dir
       basename = output_dir + options.pt +'GeV_'+str(options.proc)+'_'+str(job)
       basename = os.path.basename(basename)
       print "output name : ",basename
       currentDir = os.getcwd()
       inputFile = list_of_files[job]
              
       outputFile = output_dir+'/'+basename+'.root'

       subFile = condorPath+'/sub/'+basename+'.sub'
       exFile = condorPath+'/sub/'+basename+'.sh'
       
       fsub = None
       frun = None
       try:
           fsub = open(subFile, 'w')
       except IOError as e:
           print "I/O error({0}): {1}".format(e.errno, e.strerror)
           time.sleep(10)
           fsub = open(subFile, 'w')
       print fsub
       frun = None

       try:
           frun = open(exFile, 'w')
       except IOError as e:
           print "I/O error({0}): {1}".format(e.errno, e.strerror)
           time.sleep(10)
           frun = open(exFile, 'w')
       print frun
       commands.getstatusoutput('chmod 777 %s'%(subFile))
       
       frun.write('#!/bin/bash\n')
       frun.write(currentDir+'/submitJets.sh '+inputFile+' '+outputFile+' '+str(max_events)+' '+deltaR+' '+etaCut+' '+str(check)+' \n')
       frun.close()

       commands.getstatusoutput('chmod 777 %s'%(exFile))
       commands.getstatusoutput('chmod 777 %s'%(condorPath+'/sub/'+basename+'.sh'))
       
       fsub.write('executable            = %s\n'%(exFile))
       fsub.write('arguments             = $(ClusterID) $(ProcId)\n')
       fsub.write('output                = %s/out/job.%s.$(ClusterID).$(ProcId).out\n'%(condorPath,basename))
       fsub.write('log                   = %s/log/job.%s.$(ClusterID).$(ProcId).log\n'%(condorPath,basename))
       fsub.write('error                 = %s/err/job.%s.$(ClusterID).$(ProcId).err\n'%(condorPath,basename))
       fsub.write('RequestCpus = 4\n')
       fsub.write('+JobFlavour = "nextweek"\n')
       fsub.write('+AccountingGroup = "group_u_FCC.local_gen"\n')
       fsub.write('queue 1\n')
       fsub.close()
       
       cmd = "condor_submit %s\n"%(subFile)
      
       print cmd
       batchid=-1
       jobb,batchid=SubmitToCondor(cmd,10)       # submitting jobs
       
#       os.system(cmd)


#_______________________________________________________________________________________
if __name__ == "__main__":
    main()
