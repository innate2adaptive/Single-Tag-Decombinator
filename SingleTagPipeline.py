import os
import sys
import collections

from SingleTagDecombinator import args, get_chain
#import SupplementaryScripts.SingleTagTools.reconstructTCR as reconstructTCR
import reconstructTCR
from argparse import Namespace

import reconstructTCR

def getTagFolder():
	import urllib2
	try:
		d = "https://raw.githubusercontent.com/innate2adaptive/Decombinator-Tags-FASTAs/master/"
		urllib2.urlopen(urllib2.Request(d))      # Request URL, see whether is found
		print "here3"
	except:
		cwd = os.getcwd()
		basedir = os.path.dirname(cwd)
		if os.path.isdir(cwd+os.sep+"Decombinator-Tags-FASTAs"):
			return cwd+os.sep+"Decombinator-Tags-FASTAs"
		elif os.path.isdir(basedir+os.sep+"Decombinator-Tags-FASTAs"):
			return basedir+os.sep+"Decombinator-Tags-FASTAs"
		else:
			print "Error: Cannot find online or offline version of Decombinator-Tags-FASTAs directory."
			sys.exit()

def organiseOutput():

	dirname = "SingleTagAnalysis"
	dirnum = ""
	counter = 0
	exitflag = False
	while not exitflag:
		dirout = dirname+dirnum
		if not os.path.exists(dirout):
			os.makedirs(dirout)
			exitflag = True
		counter += 1
		dirnum = str(counter)
	return(dirout)
	

def pprint(text):
	stext = text.split(" ")
	for i in range(len(stext)):
		if i%5 == 0:
			stext.insert(i,"\n")
	print " ".join(stext)+"\n"
	return 0

def getOutputFile(inputargs):
	# Code copied from SingleTagDecombinator to obtain correct output file name

	chain = args.chain.split(" ")
	chainnams = {"a": "alpha", "b": "beta", "g": "gamma", "d": "delta"}
	inner_filename_chains = [x for x in chainnams.values() if x in inputargs.fastq.lower()]
	chains_detected = len(inner_filename_chains)
	suffix = "." + inputargs.extension

	bnam1 = os.path.basename(inputargs.fastq)
	snam1 = bnam1.split(".")[0]
	if inputargs.fastq2: #naming hack if two reads are included in input. Works with currently name data files separated by "_"
		bnam2 = os.path.basename(inputargs.fastq2)
		snam2 = bnam2.split(".")[0]
		samplenam = "_".join(collections.OrderedDict.fromkeys((snam1+"_"+snam2).split("_")).keys())
	else:
		samplenam = snam1

	if chains_detected == 1:
		name_results = inputargs.prefix + samplenam
	else:
		name_results = inputargs.prefix + "_".join(map(chainnams.__getitem__, chain)) + "_" + samplenam
	if inputargs.dontgzip == False:
		outfilenam = name_results + suffix + ".gz"
	else:
		outfilenam = name_results + suffix
	
	return outfilenam

def getDcrScript():
	import urllib2
	try:
		f = "https://raw.githubusercontent.com/innate2adaptive/Decombinator/93614c643b9151a6901f263bbb7894e46782ce07/Decombinator.py"
		urllib2.urlopen(urllib2.Request(f))      # Request URL, see whether is found
		return "curl "+f+" | python - "
		print ""
	except:
		cwd = os.getcwd()
		basedir = os.path.dirname(cwd)
		if os.path.isfile(cwd+os.sep+"Decombinator.py"):
			return "python "+cwd+os.sep+"Decombinator.py "
		elif os.path.isfile(basedir+os.sep+"Decombinator/Decombinator.py"):
			return "python "+basedir+os.sep+"Decombinator/Decombinator.py "
		elif os.path.isfile(basedir+os.sep+"Decombinator.py"):
			return "python "+basedir+os.sep+"Decombinator.py "
		else:
			print "Error: Cannot find online or offline version of Decombinator."
			sys.exit()


def Decombinator(dcr_args,outputfiles):
	cmd = getDcrScript()

	for c in dcr_args.chain.split(" "):
		dcr_input = cmd
		for a in vars(dcr_args):
			if a == 'chain':
				dcr_input += " "+"--"+"chain "+c
			elif vars(dcr_args)[a] == True:
				dcr_input += " "+"--"+a
			elif vars(dcr_args)[a] != None and vars(dcr_args)[a] != False:
				dcr_input += " "+"--"+a+" "+"\'"+str(vars(dcr_args)[a])+"\'"

		print "\n###############################"
		print "Running Decombinator: chain", c
		print "###############################\n"
		pprint("Command Issued: "+dcr_input)
		os.system(dcr_input)
		outname = getOutputFile(dcr_args)
		os.rename(outname, outdir+os.sep+c+"_"+outname)
		outputfiles.append(outdir+os.sep+c+"_"+outname)

		if dcr_args.nobarcoding == True:
			nbcfile = os.path.splitext(outname)[0]+".nbc"
			os.rename(nbcfile, outdir+os.sep+c+"_"+nbcfile)
			outputfiles.append(outdir+os.sep+c+"_"+nbcfile)


if __name__ == '__main__':

	args = args()	
	software_dir = os.path.dirname(__file__)
	#args.tagfastadir = getTagFolder()

	outdir = organiseOutput()
	outputfiles = []

	st_dcr_input = "python " + software_dir + "/SingleTagDecombinator.py"	

	for a in vars(args):
		if vars(args)[a] == True:
			st_dcr_input += " "+"--"+a
		elif vars(args)[a] != None and vars(args)[a] != False:
			st_dcr_input += " "+"--"+a+" "+"\'"+str(vars(args)[a])+"\'"


	print "###############################"
	print "Running Single Tag Decombinator"
	print "###############################"
	pprint("Command Issued: "+st_dcr_input)
	
	os.system(st_dcr_input)
	outname = getOutputFile(args)
	os.rename(outname, outdir+os.sep+outname)
	outputfiles.append(outdir+os.sep+outname)


	# recon_args = Namespace(buildfordecombinator = True,
	# 					   chains = ['a', 'b', 'c', 'd'],
	# 					   filename = outdir+os.sep+outname,
 #                           minoverlap = 5,
 #                           nocollapse = False,
 #                           noreconstruct = False,
 #                           outputfile = outdir+os.sep+"OUTPUT.n12",
 #                           separatedir = None)

 	recon_args = Namespace(filename = outdir+os.sep+outname, nproc = None)

	print "\n################################################"
	print "Running ReconstructTCR: Building For Decombinator"
	print "################################################\n"

	# bfdname = reconstructTCR.main(recon_args)
	bfdname = reconstructTCR.main(recon_args)
	os.rename(bfdname, outdir+os.sep+bfdname)
	outputfiles.append(outdir+os.sep+bfdname)
#	outputfiles.append(recon_args.outputfile)
	
	dcr_args = Namespace(allowNs = args.allowNs, 
						 bclength = 42, 
						 chain = args.chain, 
						 dontcheck = args.dontcheck, 
						 dontcount = args.dontcount, 
						 dontgzip = args.dontgzip, 
						 extension = 'n12', 
						 fastq = outdir+os.sep+bfdname, 
						 fastq2 = None,
					     lenthreshold = args.lenthreshold, 
					     nobarcoding = args.nobarcoding, 
					     orientation = 'both', 
					     prefix = 'dcr_', 
					     species = args.species, 
					     suppresssummary = args.suppresssummary, 
					     tagfastadir = None,
					     tags = args.tags)

	Decombinator(dcr_args,outputfiles)

	print "\n#######################################################"
	print "Output Files have been saved to:"
	for f in outputfiles:
		print f
	print "########################################################\n"