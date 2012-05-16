import re
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio import SeqIO
from optparse import OptionParser


#Parse command line options
parser = OptionParser()
parser.add_option("-w", "--input", dest="wellId", action="store", help=".Well ID", metavar="STRING")
parser.add_option("-o", "--output", dest="outputTag", help="output file tag", metavar="STRING")
(options,args) = parser.parse_args()
revReadRegEx = re.compile('.*R.fq')
primerIds = ["2_4_F.fq","2_4_R.fq","2_2_F.fq","2_5_R.fq",]
wellId="/data/A04/MiSeq/%s/"%options.wellId
dnaRegEx = re.compile('^[ACGTN]+$')



def reverse_complement(read):
    seq = Seq(read, generic_dna)
    return seq.reverse_complement().tostring()    

def rank2(finalCounts, recur):
    readTotals = [0]*4
    counter=0
    for entry in finalCounts.items():
        for x in range (0,4):
            readTotals[x]+=entry[1][x]
    
    
    countsPerRead = [0] *4
    for readNum in range(0,4):
        #print readNum
        countsPerRead[readNum]={}
        for entry in finalCounts.items():
            name = entry[0]
            counts = entry[1]
            countsPerRead[readNum][name] = counts[readNum]
            #print readNum,countsPerRead[readNum][name]
        
    #now rank each gene based on number of reads mapping
    finalRanks = {}
    for readNum in range(0,4):
        #print "Analysing readNum %d"%readNum
        #read through all readCounts in sorted order and assign a rank
        currRank =1
        prevScore = -1
        allSameScoreList = []
        rankDifferential = 0 #this is so that if you had scores of 10, 20, 20, 50, t
        for entry in sorted(countsPerRead[readNum], key=countsPerRead[readNum].get, reverse=True): 
            score = countsPerRead[readNum][entry]
            name = entry 
            if score != prevScore:
                currRank = currRank+rankDifferential+1      
                rankDifferential=1
                allSameScoreList = []
                allSameScoreList.append(name)
                #add rank
                if name not in finalRanks:
                    finalRanks[name] = [0]*4
                finalRanks[name][readNum]= currRank                
            else:
                #add rank
                rankDifferential+=1
                if name not in finalRanks:
                    finalRanks[name] = [0]*4
                finalRanks[name][readNum]= rankDifferential + currRank  
                allSameScoreList.append(name)
                for name in allSameScoreList:
                    finalRanks[name][readNum]=rankDifferential + currRank
            prevScore = score
            #rankDict[readNum][name]=currRank
#            if name not in finalRanks:
#                finalRanks[name] = [0]*4
#            finalRanks[name][readNum]= currRank
#           optimisticRank[name
           # print name,score,readNum,currRank, finalRanks[name],finalCounts[name],[readNum]
            
    #calculate total rankings
    overallRanks = {}
    for entry in finalRanks.items():
        #print entry
        overallRanks[entry[0]] = entry[1][0] + entry[1][1] + entry[1][2] + entry[1][3]
    
    newScoreMap = {}
    for entry in sorted(overallRanks, key=overallRanks.get, reverse=True):
        name = entry
        rank = overallRanks[entry]
        individualRanks = finalRanks[name]
        originalCounts = finalCounts[name]

        
        #new ranking idea multiply highest count by worst rank.....lowest count by best rank
        orderedRanks=[]
        orderedCounts=[]
        for i in sorted(individualRanks):
            orderedRanks.append(i)
        for i in sorted(originalCounts):
            orderedCounts.append(i)
        newIdeaScore = 0.0
        for i in range(0,4):
            newIdeaScore+=orderedCounts[i]/orderedRanks[i]
        
        misses =0
        calibratedScore = 0
        missModifier = 1.0        
        for i in range(0,4):
            calibratedScore+=originalCounts[i]/individualRanks[i]
            if originalCounts[i]==0:
                missModifier-=0.25
                misses+=1
                
        newIdeaScore = newIdeaScore/(1+misses)
        calibratedScore = calibratedScore / missModifier
        newScoreMap[name]=newIdeaScore
    #    print "size",len(newScoreMap), name, newIdeaScore, newScoreMap[name]
        #add to map based on newIdeaScore
        
        print name,"\t",individualRanks,"\t",originalCounts,"\t", rank,"\t", sum(originalCounts)/rank,"\t", calibratedScore,"\t", newIdeaScore
    
    top10 =0
   # for entry in newScoreMap.items():
      #  print "@",entry
    #for entry in sorted(finalScores, key=finalScores.get, reverse=True):
    for entry in sorted(newScoreMap, key=newScoreMap.get, reverse=True):
      #  print "ENTRY:",entry, len(newScoreMap)
        name = entry
        score = newScoreMap[entry]
        
        if top10 < 10:
            print name,finalRanks[name],finalCounts[name],score
        top10+=1
            
    
    #translate back to score map i.e. (name,[1][10][13][1])
#    for entry in finalRanks.items():
#        print entry[0],entry[1]
            


def rank(finalCounts, recur):
    readTotals = [0]*4
    counter=0
    for entry in finalCounts.items():
#        print entry
        for x in range (0,4):
            readTotals[x]+=entry[1][x]
    
#    print "ReadTotals:%s" % readTotals
    
    overallTotals = sum(readTotals)
    readProportions = [0.0]*4
    for x in range(0,4):
        readProportions[x] = readTotals[x]/float(overallTotals)
#    print "ReadProps :%s" % readProportions
    
    #now, calculate all scores
    finalScores = {}
    for entry in finalCounts.items():
        noHits=0
        name = entry[0]
        counts = entry[1]
        normalisedCounts = [0.0]*4
        totalForTarget = sum(counts)
        for x in range(0,4):
            if counts[x]==0:
        #        print "NOHITS:",name,noHits
                noHits+=1
            if totalForTarget!=0 and readProportions[x]>0.0001:
                normalisedCounts[x] = counts[x]/(totalForTarget*readProportions[x])
         #   print x,normalisedCounts[x]
            if normalisedCounts[x] > 1.0:
                normalisedCounts[x] = 1.0
       
        multipliedScore = (normalisedCounts[0]+0.01) * (normalisedCounts[1]+0.01) * (normalisedCounts[2]+0.01) * (normalisedCounts[3]+0.01)
        finalScore = 0
        if totalForTarget!=0:
            finalScore = pow(multipliedScore,(1.0/float(totalForTarget)))
        
        #normalise based on missing read mappings
        readMissingNormaliser = 1 / (1+ pow(noHits,0.1)/1000)
        
        finalScores[name]=finalScore * readMissingNormaliser# float(1.0/(1.0+((noHits*noHits)/100.0))) #float(1.0/(1+(noHits/10.0)))
        #TODO: This may elimiate otherwise good candidates - therefore MUST have a warning system in place - if the score is high but gets demoted due to a read not mapping to it, MUST report it's pre-calibration score too!
       # print name, finalScore, finalScore* float(1.0/(1.0+((noHits*noHits)/100.0)))
    
   
    #only bother printing out if this is the last recurrance of rank
    if recur==False:
        top10=0
        print
        print "SCORES",recur
        
        for entry in sorted(finalScores, key=finalScores.get, reverse=True):
                score = finalScores[entry]
                name = entry
                #tallies = counts
                tallies = finalCounts[name]
                if top10 < 10:
                    print name,tallies, score  
                top10+=1                
  
    #rebuild a suitable construct to send back to this function - this will recalibrate based only on the top 10 results
    RECAL_BASED_ON = 30
    if recur == True:
        top10=0
        recalMap = {}
        for entry in sorted(finalScores, key=finalScores.get, reverse=True):    
            score = finalScores[entry]
            name = entry
            tallies = finalCounts[name]
            if top10 < RECAL_BASED_ON:
                recalMap[name]=tallies
            top10+=1
        
        rank2(recalMap,False)
        rank(recalMap,False)
        
    
    
    #return top 10
    top10=0
    returnMap = {}
    for entry in sorted(finalScores, key=finalScores.get, reverse=True):    
        score = finalScores[entry]
        name = entry
        if top10 < 300000:
            returnMap[name]=score
        top10+=1    
    
    return returnMap
        
                         
                        
            

countsDict = {}
countsDict["2_2_F.fq"]={}
countsDict["2_4_F.fq"]={}
countsDict["2_4_R.fq"]={}
countsDict["2_5_R.fq"]={}

#load all ref seqs into memory
refs = {}
antiRefs = {}
drb1RegEx = re.compile('DRB1')
refFile = open("subject.txt")
for line in refFile:
    line = line.strip()
    
    lineSpl = line.split("\t")
    if len(lineSpl) ==2:
        if drb1RegEx.match(line):
            refs[lineSpl[0]] = lineSpl[1]
        else:
            antiRefs[lineSpl[0]] = lineSpl[1]   

refFile.close()

for p in primerIds:
    currCounts = {}
    readCount =0
    hits =0
    #change behaviour based on which read it is
    if p== "2_2_F.fq": #0:160
        reverse = False
        startPos = 0
        endPos = 160
    elif p== "2_4_F.fq":#0:150
        startPos = 0
        endPos = 150
        reverse = False
    elif p=="2_4_R.fq":#50:120
        startPos = 33
        endPos = 200#250
        reverse = False
    elif p=="2_5_R.fq": #70:230
        startPos = 70
        endPos = 232
        reverse = False   
       
       
    queryFile = open("%s%s"%(wellId,p))
    outFile = open("%s_results.txt"%(p),"w")
    for line in queryFile:
        
        if dnaRegEx.match(line):  
            query= line[startPos:endPos]
            if reverse == True:
                query = reverse_complement(query)
           
            readCount+=1
            presentInAntiRef = False
            for aRef in antiRefs.items():
                aName = aRef[0]
                aSeq = aRef[1]
                if query in aSeq:
                    #print "Found in %s"%aName
                    presentInAntiRef = True
                    
                    for ref in refs.items():
                        name = ref[0]
                        seq = ref[1]
                      #  print name,seq
                    #    if query in seq:   
                    #        print "Found in %s and %s"%(aName,name)
                    break
            if presentInAntiRef == False:
                #print query
                for ref in refs.items():
                    name = ref[0]
                    seq = ref[1]
                   # print name,seq
                    if query in seq:
                      #  print "\t",name
                        if name not in currCounts:
                            currCounts[name] = 1
                        else:
                            currCounts[name] +=1
                        hits+=1

    outFile.close()
    queryFile.close()
    
    countsDict[p]=currCounts
    
#print currCounts


#now create counts based on only first 2 fields of name
for p in primerIds:
    shorterCounts = {}
    for hit in countsDict[p].items():
        colonSpl = hit[0].split(":")
        shorterName = "%s:%s"%(colonSpl[0],colonSpl[1])    
        if shorterName not in shorterCounts:
            shorterCounts[shorterName]=hit[1]
        else:
            shorterCounts[shorterName]+=hit[1]
    countsDict[p]=shorterCounts
        
#compare counts
#add all detected hlas to a set
allHits = set()
for p in primerIds:
    for hit in countsDict[p].items():
        allHits.add(hit[0])
        
#print "ALL HITS: %s"%allHits

# now do counts per hit across samples present in all 

#print
#print
finalMap = {}
success = False
for commonHit in allHits:
    
    finalMap[commonHit]=[0]*4
    currSample = 0
    for p in primerIds:
        if commonHit in countsDict[p]:
            finalMap[commonHit][currSample] = countsDict[p][commonHit]
          #  print commonHit, countsDict[p][commonHit]
        currSample+=1
#    
#    
#    
#    if 0 not in finalMap[commonHit]:
#        success = True
#        print "%s\t%d\t%d\t%d\t%d\t%d"%(commonHit,finalMap[commonHit][0],finalMap[commonHit][1],finalMap[commonHit][2],finalMap[commonHit][3],sum(finalMap[commonHit]))

#if success == False:
#    for entry in finalMap.items():
#        commonHit = entry[0]
#        zeroCount=0
#        for c in finalMap[commonHit]:
#            if c==0:
#                zeroCount+=1
#        if zeroCount <=2:
#            print "%s\t%d\t%d\t%d\t%d\t%d"%(commonHit,finalMap[commonHit][0],finalMap[commonHit][1],finalMap[commonHit][2],finalMap[commonHit][3],sum(finalMap[commonHit]))


finalScores = rank(finalMap, True)

#remap(finalScores,refs)


        
    