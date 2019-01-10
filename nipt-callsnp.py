# coding=utf-8
import os,datetime,shutil,subprocess,time
import psutil#,threading
import multiprocessing 
#import threading
#from threading import Thread

os.chdir('/home/userroot/data/nipt')
files = open('./call_snp/C6R5.csv')
fq_dir = files.read().splitlines()
files = open('./call_snp/fq_C6R5.csv')
fq_name = files.read().splitlines()


i = datetime.datetime.now()
start_time = str(i)[0:19]
#print(fq_name)
filename_fq=open('./filename/'+str(i.year)+'-'+str(i.month)+'-'+str(i.day)+'-'+str(i.hour)+'：'+str(i.minute)+'.txt',"w")
log=open('./log/'+str(i.year)+'-'+str(i.month)+'-'+str(i.day)+'-'+str(i.hour)+'：'+str(i.minute)+'.txt',"a")
filename_fq.writelines([line+'\r\n' for line in fq_name])
n = len(fq_name)
m = n

#######>> making bin file and gc_bin file <<########
#os.system("bedtools makewindows -g ./bed/chr_size/chr_size.sizes -w 300000 > ./bed/bin_ucsc_hg19.bed")
#os.system("bedtools nuc -fi ./fa/fasta/ucsc.hg19.fasta -bed ./bed/bin_ucsc_hg19.bed | cut -f 1-3,5 > ./gc_bed/gc_bin.bed")

class nipt(object):
	'''parameter_1: fq files'''
	'''parameter_2: direction of fq files(keyboard input)'''
	'''cpu rate cutoff and wait time is set for 96 or more samples.DON'T CHANGE!'''
	def __init__(self,fq,file_dir):
		self.fq = fq
		self.file_dir = file_dir
		

	def setdir(self):
		#input_file= './fq/'+self.file_dir+'/'+self.fq
		input_file= '/home/cyto/nipt_back_up/fq/'+self.file_dir
		sam_file= './sam/'+self.fq+'.sam'
		bam_file= './bam/'+self.fq+'.bam'
		sort_bam_file='./sort_bam/sort_'+self.fq+'.bam'
		#sort_bam_file_1='./sort_bam/sort_'+self.fq+'.bam'
		rmdup_sort_bam_file= 'rmdup_sort_bam/rmdup_sort_'+self.fq+'.bam'
		#counts_file='./counts/'+self.fq+'.counts' 
		gvcf_file = './gvcf/'+self.fq+'.g.vcf'
		file_id = self.fq
		METRICS = 'rmdup_sort_bam/'+self.fq+'.metrics'
		PL = 'illumina'
		LB = 'NIPT'
		SM = '01'
		a = [input_file,sam_file,bam_file,sort_bam_file,rmdup_sort_bam_file,gvcf_file,file_id,PL,LB,SM,METRICS]
		return a
		
	def mapping(self,zuse = 0):
		a = self.setdir()
		print(a[0])
		print(a[1])
		root_dir = '/home/userroot/data/nipt/'
		RG = '@RG\\tID:'+a[6]+'\\tSM:'+a[6]+'\\tLB:WES\\tPL:ILLUMINA'
		print(RG)
		b = 'bwa mem -t 40 -R "'+RG+'" -o ' + a[1] + ' ./fa/ucsc_hg38_bwa/Homo_sapiens_assembly38.fasta ' +a[0]
		c = 'bwa mem -t 40  -o ' + a[1] + ' ./fa/ucsc_hg38_bwa/Homo_sapiens_assembly38.fasta ' +a[0]
		#print(b)
		p = subprocess.Popen([b],shell=True)
		#p = subprocess.Popen(["subread-align -t 1 -T 40 -M 0 -I 0 -i ./fa/ucsc_hg19/ucsc_hg19 -r"+a[0]+" -o"+a[1]],shell=True) #mapping
		if zuse == 1:
			p.wait()
		return 
		
	def samtobam(self):
		a = self.setdir()
		print(a[1])
		print(a[2])
		#subprocess.Popen("rm ./bam/*.indel",shell=True) 
		p = subprocess.Popen(["samtools view -b -@ 40 "+a[1]+" > "+a[2]],shell=True) #sort
		p.wait()
		#os.remove(a[1])
		return  
		
	def sorting(self):
		a = self.setdir()
		#subprocess.Popen("rm ./bam/*.indel",shell=True) 
		p = subprocess.Popen(['~/software/gatk-4.0.4.0/gatk SortSam -I '+a[2]+ ' --SORT_ORDER coordinate -O ' +a[3]],shell=True)
		#p = subprocess.Popen(["samtools sort -@ 40 -o "+a[2]+" "+a[1]],shell=True) #sort
		#p2.wait()
		#os.remove(a[1])
		return  
  
	def rmdup(self):
		a = self.setdir()
		p = subprocess.Popen(['~/software/gatk-4.0.4.0/gatk MarkDuplicates -I '+ a[3] + ' --REMOVE_DUPLICATES true --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 8000 -O '+a[4] +' --METRICS_FILE '+ a[10]],shell=True)
		#p = subprocess.Popen(["samtools rmdup -s "+a[3]+" "+a[4]],shell=True) #remove duplications
		#p3.wait()
		#os.remove(a[3])
		return 
  
	def index(self):
		a = self.setdir()
		p = subprocess.Popen('samtools index -@ 40 '+a[4],shell=True) 
		p.wait()
		#os.remove(a[4])
		return
		
	def snp(self):
		a = self.setdir()
		p = subprocess.Popen(['~/software/gatk-4.0.4.0/gatk HaplotypeCaller -ERC GVCF  -R /home/userroot/data/nipt/fa/hg38/Homo_sapiens_assembly38.fasta -D /home/userroot/data/nipt/annotation/dbsnp_146.hg38.vcf  -I /home/userroot/data/nipt/'+a[4] +' -O ' +a[5]],shell=True)
		#p = subprocess.Popen('samtools index -@ 40 '+a[4],shell=True) --dbsnp /home/userroot/data/nipt/annotation/dbsnp_146.hg38.vcf
		#p4.wait()
		#os.remove(a[4])
		return   
 
	def align(self,zuse = 0):
		cpu = psutil.cpu_percent(interval=1)
		mem = psutil.virtual_memory()
		print(mem.percent)
		if float(cpu) < 65.0 and mem.percent < 30:
			self.setdir()
			self.mapping()
			time.sleep(8)
		else:
			time.sleep(5)
			cpu = psutil.cpu_percent(interval=1) 
			mem = psutil.virtual_memory()
			while float(cpu) > 65.0 or mem.percent > 30:
				time.sleep(8)
				cpu = psutil.cpu_percent(interval=1) 
				mem = psutil.virtual_memory()
			else:
				self.setdir()
				self.mapping()
				time.sleep(8)
				
#	def samtobam(self,zuse = 0):
#		#cpu = psutil.cpu_percent(interval=1)
#		mem = psutil.virtual_memory()
#		print(mem.percent)
#		if float(mem.percent < 80):
#			self.setdir()
#			self.sam2bam()
#			time.sleep(3)
#		else:
#			time.sleep(3)
#			#cpu = psutil.cpu_percent(interval=1) 
#			mem = psutil.virtual_memory()
#			while float( mem.percent > 80):
#				time.sleep(3)
#				#cpu = psutil.cpu_percent(interval=1) 
#				mem = psutil.virtual_memory()
#			else:
#				self.setdir()
#				self.sam2bam()
#				time.sleep(3)
						
	def sort(self):
		cpu = psutil.cpu_percent(interval=1)
		if float(cpu) < 80.0:
			self.sorting()
			time.sleep(5)
		else:
			time.sleep(5)
			cpu = psutil.cpu_percent(interval=1) 
			while float(cpu) > 80.0:
				time.sleep(5)
				cpu = psutil.cpu_percent(interval=1) 
			else:
				self.sorting()
				time.sleep(5)
				
	def callsnp(self):
		cpu = psutil.cpu_percent(interval=1)
		mem = psutil.virtual_memory()
		if float(cpu) < 80.0 and mem.percent < 60:
			self.snp()
			time.sleep(5)
		else:
			time.sleep(5)
			cpu = psutil.cpu_percent(interval=1) 
			mem = psutil.virtual_memory()
			while float(cpu) > 80.0 or mem.percent > 60:
				time.sleep(5)
				cpu = psutil.cpu_percent(interval=1) 
				mem = psutil.virtual_memory()
			else:
				self.snp()
				time.sleep(5)
				
def merge_bam():
		sort_bam = os.listdir('./rmdup_sort_bam')
		sort_bam_file = [x for x in sort_bam if 'bam' in x and 'bai' not in x]
		#print(sort_bam_file)
		sort_bam_name=list()
		for i in sort_bam_file:
			j = './rmdup_sort_bam/'+i
			sort_bam_name.append(j)
		filename_bam = open('./filename/filename_bam.txt',"w")
		filename_bam.writelines([line+'\r\n' for line in sort_bam_name])
		os.system('samtools merge -c -@ 40 -b ./filename/filename_bam.txt ./merged/C6R5.bam -O BAM')
		os.system('samtools index ./merged/C6R5.bam')
		
	
			
def onebyone(fq,filedir,method,waittime_1,waittime_2):
	for i in range(n):
		sample = nipt(fq[i],filedir[i])	
		if method == 0:
			sample.align()
		elif method == 1:
			sample.samtobam()
		elif method == 2:
			sample.sorting()
		elif method == 3:
			sample.rmdup()
		elif method == 4:
			sample.index()
		elif method == 5:
			sample.callsnp()
		#elif method == 6:
			#sample.merge_bam()
		time.sleep(waittime_1)
	time.sleep(waittime_2)
	
#########>>> This is main process <<<###########
#n = 96				 
def main():
	print(range(n))
	i = datetime.datetime.now()
	map_time = str(i)[0:19]
	onebyone(fq_name,fq_dir,0,30,300)
	i = datetime.datetime.now()
	sort_time = str(i)[0:19]
	onebyone(fq_name,fq_dir,1,0,600)
	i = datetime.datetime.now()
	rmdup_time = str(i)[0:19]
	onebyone(fq_name,fq_dir,2,45,600)
	i = datetime.datetime.now()
	counts_time = str(i)[0:19]
	onebyone(fq_name,fq_dir,3,60,600)
	i = datetime.datetime.now()
	onebyone(fq_name,fq_dir,4,0,10)
	i = datetime.datetime.now()
	merge_bam()
	i = datetime.datetime.now()
	finish_time = str(i)[0:19]
	log.write('start time: '+start_time+"\r\n"+ \
	'mapping begin at: '+map_time+"\r\n"+ \
	'sorting begin at: '+sort_time+"\r\n"+ \
	'removing dup begin at: '+rmdup_time+"\r\n"+ \
	'counting begin at: '+counts_time+"\r\n"+ \
	'finish time: '+finish_time+' \r ')
################################################

#################>>> START <<<##################
if __name__ == '__main__':
	main() 
 #for line in fq_name: 
  #t1 = threading.Thread(mapping(line))
  #threads.append(t1)
  #for t in threads:
  #  t.setDaemon(True)
  #  t.start()	
  #p = multiprocessing.Pool(processes=40)
  #for i in range(10):
  #  p.apply_async(main) 
    
    #p1 = multiprocessing.Process(mapping(line))
    #p1.start()  
    #p.apply_async(mapping(line))
  #print("#######unanalysed/total:"+str(n-1)+"/"+str(m)+";Progress status:"+str(format(1-((n-1)/m),'.0%'))+"########") #this is a counter.
  #n = n-1
 # p.close() # 开始执行进程
 # p.join() # 等待子进程结束
filename_fq.close()
log.close()



