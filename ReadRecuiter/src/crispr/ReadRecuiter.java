package crispr;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.logging.ConsoleHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;

public class ReadRecuiter {
	@Parameter(names = "-repeats", description = "Known repeats fasta file", required = false)
	private String repeatFile;

	@Parameter(names = "-minRepeat", description = "Minimum length of a DR, default=20", required = false)
	private int minRepeatLength = 20;

	@Parameter(names = "-maxRepeat", description = "Maximum length of a DR, default = 50", required = false)
	private int maxRepeatLength = 50;

	@Parameter(names = "-minSpacer", description = "Minimum length of a spacer, default = 20", required = false)
	private int minSpacerLength = 20;

	@Parameter(names = "-maxSpacer", description = "maximum length of a spacer, default = 60", required = false)
	private int maxSpacerLength = 60;

	@Parameter(names = "-maxMismatch", description = "maximum mismatch allowed in seed kmer, default = 1", required = false)
	private int maxMismatch = 1;

	@Parameter(names = "-input", description = "Input paired end reads fasta file", required = true)
	private String fastaFile;

	@Parameter(names = "-kmerSize", description = "Kmer size, default = 20", required = false)
	private int kmerSize = 20;

	@Parameter(names = "-threads", description = "Number of threads", required = false)
	private int nThread = 1;

	@Parameter(names = "--help", help = true)
	private boolean help = false;


	@Parameter(names = "-prefix", description = "Prefix for output files. Default='out", required = false)
	private String outPrefix = "out";

	private static HashSet<String> kmerSet1 = new HashSet<String>();
	private static HashSet<String> kmerSet2 = new HashSet<String>();
	private static HashMap<String, String> DRReads = new HashMap<String, String> (); 
	private static Logger logger = Logger.getLogger(ReadRecuiter.class.getName());
	public static void main(String args[]){
		ReadRecuiter settings = new ReadRecuiter();
		JCommander jc = new JCommander(settings, args);
		if (settings.help) {
			jc.usage();
			return;
		}
		//System.out.println(settings.minRepeatLength);
		//System.out.println(settings.maxRepeatLength);
		//System.out.println(settings.minSpacerLength);
		//System.out.println(settings.maxSpacerLength);

		String drFile = null;
		if (settings.repeatFile != null) {
			drFile = settings.repeatFile;
		}
		String readsFile = settings.fastaFile;
		String outPrefix = settings.outPrefix;
		int kmerSize = settings.kmerSize;
		if (settings.maxMismatch>2) {
			settings.maxMismatch = 2;
		}
		else if (settings.maxMismatch <0) {
			settings.maxMismatch = 0;
		}
		if (settings.maxMismatch != 0) {
			settings.kmerSize = settings.minRepeatLength;
		}
		int nThread = settings.nThread;
		if (nThread < 1) {
			logger.log(Level.WARNING, "Number of threads sets to 1.");
			nThread = 1;
		}

		logger.log(Level.INFO, "************************************************************.");
		logger.log(Level.INFO, "Gathering reads contain at least two copies of repeats.");

		//Creating shared object
		BlockingQueue <List<PairedEndReads>> sharedQueue = new LinkedBlockingQueue  <List<PairedEndReads>> (nThread+2);
		//Creating Producer and Consumer Thread
		Thread fastaReaderThread = new Thread(new FastaReader(sharedQueue, readsFile));
		fastaReaderThread.start();
		ArrayList<Thread> DRFinderThreads = new ArrayList<Thread> ();

		//ConcurrentHashMap<String, String> mapDR = new ConcurrentHashMap<String, String>();
		//ConcurrentHashMap<String, String> mapSeq = new ConcurrentHashMap<String, String>();
		ArrayList <String> drfiles = new ArrayList<String> ();
		for (int i=0; i<nThread; i++) {
			String outNameRepeat = outPrefix + "." + Integer.toString(i) + ".dr.fa";
			String outNameReads = outPrefix + "." + Integer.toString(i) + ".dr.reads.fa";
			drfiles.add(outNameRepeat);
			Thread t;
			try {
				t = new Thread(new RepeatReadsRecuiter(sharedQueue, settings.minRepeatLength, settings.maxRepeatLength,
						settings.minSpacerLength, settings.maxSpacerLength, settings.maxMismatch, outNameRepeat, outNameReads));
				DRFinderThreads.add(t);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} 
		}
		for (Thread t: DRFinderThreads) {
			t.start();
		}
		for (Thread t: DRFinderThreads) {
			try {
				t.join();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		logger.log(Level.INFO, "************************** DONE ************************.\n\n");


		logger.log(Level.INFO, "************************************************************.");
		logger.log(Level.INFO, "Gathering reads with one copy of repeat.");

		logger.log(Level.INFO, "Generating kmers from recuited repeats.");
		for (String name : drfiles) {
			genKmersFromDBFile(name, kmerSize);
		}
		logger.log(Level.INFO, String.format("Kmer hash table size: %d", kmerSet1.size()));

		if (drFile != null) {
			logger.log(Level.INFO, "Generating kmers from supplied known repeats.");
			genKmersFromDBFile(drFile, kmerSize);
			logger.log(Level.INFO, String.format("Kmer hash table size: %d", kmerSet1.size()));
		}

		//Creating shared object
		sharedQueue = new LinkedBlockingQueue  <List<PairedEndReads>> (nThread+2);
		//Creating Producer and Consumer Thread
		fastaReaderThread = new Thread(new FastaReader(sharedQueue, readsFile));
		fastaReaderThread.start();
		ArrayList<Thread> filterThreads = new ArrayList<Thread> ();
		for (int i=0; i<nThread; i++) {
			String outName = outPrefix + "." + Integer.toString(i) + ".filtered.fa";
			String outNameNameMapping = outPrefix + "." + Integer.toString(i) + ".filtered.name.des";
			Thread t;
			try {
				t = new Thread(new ReadFilter(sharedQueue, kmerSize, ReadRecuiter.kmerSet1, ReadRecuiter.kmerSet2, DRReads, outName, outNameNameMapping, i, nThread));
				filterThreads.add(t);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		for (Thread t: filterThreads) {
			t.start();
		}

		logger.log(Level.INFO, "************************** DONE ************************.\n\n");

	}

	private static List <String> getKmerFromString(String seq, int kmerSize) {
		ArrayList<String> kmers = new ArrayList<String> ();
		if (kmerSize > seq.length()) 
			return kmers;
		int end = seq.length() - kmerSize;
		for (int i=0; i<=end; i++) {
			kmers.add(seq.substring(i, i+kmerSize));
		}
		return kmers;
	}

	public static char reverseCompliment(char c) {
		if (c == 'A')
			return 'T';
		if (c == 'G')
			return 'C';
		if (c == 'T')
			return 'A';
		else
			return 'G';
	}

	public static String reverseCompliment(String seq) {
		StringBuilder sb = new StringBuilder();
		for (int i=seq.length()-1; i>=0; i--) {
			sb.append(reverseCompliment(seq.charAt(i)));
		}
		return sb.toString();
	}

	private static void genKmersFromDBFile(String inputfile, int kmerSize) {
		String line;
		try {
			File file = new File(inputfile);
			FileReader fileReader = new FileReader(file);
			BufferedReader bufferedReader = new BufferedReader(fileReader);
			while ((line = bufferedReader.readLine()) != null) {
				if (!line.startsWith(">")) {
					List<String> kmers = getKmerFromString(line, kmerSize);
					for (String e: kmers) {
						kmerSet1.add(e);
					}
					kmers = getKmerFromString(reverseCompliment(line), kmerSize);
					for (String e: kmers) {
						kmerSet2.add(e);
					}
				}
				else {
					String [] sp = line.split("\\t");
					DRReads.put(sp[0], sp[1]);
				}
			}
			bufferedReader.close();
		}catch (Exception e) {
			e.printStackTrace();
		}
	}
}


class FastaReader implements Runnable {
	private final BlockingQueue <List<PairedEndReads>> sharedQueue;
	private int blockSize = 50000;
	private BufferedReader bufferedReader;
	private int queueCapacity;
	private Logger logger = Logger.getLogger(FastaReader.class.getName());
	public FastaReader(BlockingQueue <List<PairedEndReads>> sharedQueue, String fastaFile) {
		this.sharedQueue = sharedQueue;
		this.queueCapacity = sharedQueue.remainingCapacity();
		File file = new File(fastaFile);
		FileReader fileReader;
		try {
			fileReader = new FileReader(file);
			this.bufferedReader = new BufferedReader(fileReader);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	@Override
	public void run() {
		try {
			List<PairedEndReads> readlist;
			while ((readlist = readBlock()) != null) {
				logger.log(Level.INFO, "Reading a block...");
				sharedQueue.put(readlist);
				logger.log(Level.INFO, "Done read a block...");
			}
			logger.log(Level.INFO, "All data read, add poison pills..");
			// poison pills
			List<PairedEndReads> pill = new ArrayList<PairedEndReads> ();
			for (int i=0; i < queueCapacity; i++) {
				logger.log(Level.INFO, "Add pill");
				sharedQueue.put(pill);
			}
		} catch (InterruptedException ex) {
			logger.log(Level.SEVERE, null, ex);
		}
	}

	private List<PairedEndReads> readBlock() {
		List<PairedEndReads> readList = new ArrayList<PairedEndReads> ();
		String id1;
		String seq1;
		String id2;
		String seq2;
		try {
			while ((id1 = bufferedReader.readLine()) != null) {
				id1 = id1.split("\\s+")[0];
				seq1 = bufferedReader.readLine();
				id2 = bufferedReader.readLine().split("\\s+")[0];
				seq2 = bufferedReader.readLine();
				PairedEndReads r = new PairedEndReads(id1, seq1, id2, seq2);
				readList.add(r);
				if (readList.size() == blockSize) {
					return readList;
				}
			}
			if (readList.size() > 0) {
				return readList;
			}
			else
				return null;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}

}


class RepeatReadsRecuiter implements Runnable{
	class RepeatInfo {
		private String repeat;
		private List<Integer> positions;
		private boolean isValid;
		public RepeatInfo() {
			positions = new ArrayList<Integer> ();
			isValid = false;
		}
		public void setRepeat(String r) {
			this.repeat = r;
		}
		public String getRepeat() {
			return this.repeat;
		}
		public void addPosition(int pos) {
			this.positions.add(pos);
		}
		public List<Integer> getPositions() {
			return this.positions;
		}
		public void setValid(boolean valid) {
			this.isValid = valid;
		}
		public boolean isValid () {
			return this.isValid;
		}
		public void reset() {
			this.repeat = null;
			this.positions.clear();
			this.isValid = false;
		}
	}

	private int minRepeatLength;
	private int maxRepeatLength;
	private int minSpacerLength;
	private int maxSpacerLength;
	private int maxMismatch;
	private BufferedWriter bufferedWriterDr;
	private BufferedWriter bufferedWriterReads;
	private final BlockingQueue <List<PairedEndReads>>  sharedQueue;
	private int seedSize = 12;
	private Logger logger = Logger.getLogger(RepeatReadsRecuiter.class.getName());
	public RepeatReadsRecuiter (BlockingQueue <List<PairedEndReads>> sharedQueue,
			int minRepeatLength, int maxRepeatLength,
			int minSpacerLength, int maxSpacerLength, int maxMismatch, String outNameDr, String outNameReads) throws IOException {
		this.sharedQueue = sharedQueue;
		this.minRepeatLength = minRepeatLength;
		this.maxRepeatLength = maxRepeatLength;
		this.minSpacerLength = minSpacerLength;
		this.maxSpacerLength = maxSpacerLength;
		this.maxMismatch = maxMismatch;
		this.bufferedWriterDr = new BufferedWriter(new FileWriter(new File(outNameDr)));
		this.bufferedWriterReads = new BufferedWriter(new FileWriter(new File(outNameReads)));
		this.logger.setLevel(Level.ALL);
		Handler consoleHandler = new ConsoleHandler();
		//consoleHandler.setLevel(Level.FINER);
		this.logger.addHandler(consoleHandler);
		if (this.maxMismatch==0) {
			seedSize = this.maxMismatch;
		}
	}

	@Override
	public void run() {
		int readsProcessed = 0;
		while(true){
			try {
				logger.log(Level.INFO, "Taking a new list");
				List<PairedEndReads> inList = sharedQueue.take();
				logger.log(Level.INFO, "New list size: " + Integer.toString(inList.size()));
				readsProcessed += 2 * inList.size();
				if (inList.size() == 0) {
					logger.log(Level.INFO, "Done. RepeatReadsRecuiter thread exit..");
					try {
						this.bufferedWriterDr.close();
						this.bufferedWriterReads.close();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					break;
				}
				else {
					int withAtLeastTwoDR = checkReadsForDR(inList);
					float percent = 100 * (float)(withAtLeastTwoDR) / inList.size();
					logger.log(Level.INFO, String.format("One block done. Total reads proccessed by this thread:  %d reads, recuited %.3f pecent of the input reads",
							readsProcessed, percent));
				}
			} catch (InterruptedException ex) {
				logger.log(Level.SEVERE, null, ex);
			}
		}
	}

	private boolean containOnlyOneChar(String s) {
		char first = s.charAt(0);
		for (int i=1; i<s.length(); i++) {
			if (s.charAt(i) != first) {
				return false;
			}
		}
		return true;
	}

	private boolean goodRepeat(String s) {
		float a = 0;
		float g = 0;
		float c = 0;
		float t = 0;
		for (int i=0; i<s.length(); i++) {
			char ch = s.charAt(i);
			switch (ch) {
			case 'A':
				a++;
				break;
			case 'G':
				g++;
				break;
			case 'C':
				c++;
				break;
			default:
				t++;
			}
		}
		a = a/s.length();
		g = g/s.length();
		c = c/s.length();
		t = t/s.length();
		if (a > 0.85 ||g > 0.85 || c > 0.85 ||t > 0.85) {
			return false;
		}
		return true;
	}

	private ArrayList<Integer>  allIndex(String seq, String pattern) {
		ArrayList<Integer> ret = new ArrayList<Integer> ();
		int lastIndex = -1;
		while (true) {
			lastIndex = seq.indexOf(pattern, lastIndex+1);
			if (lastIndex>=0) {
				ret.add(lastIndex);
			}
			else {
				break;
			}
		}
		return ret;
	}

	private boolean goodSeq(String seq, int p, int q, int len) {
		// the sequence from the first repeat to the second repeat should only contain
		// two copies of the repeat sequence. 		
		String repeat1 = seq.substring(p, len+p);
		String repeat2 = seq.substring(q, len+q);
		if (repeat1.compareTo(repeat2)==0) {
			logger.log(Level.FINE, "No mismatch in repeats");
			ArrayList<Integer> index = allIndex(seq.substring(p, q+len), repeat1);
			if (index.size() == 2) {
				return true;
			}
			else {
				logger.log(Level.FINE, "More than two repeats:" + Integer.toString(index.size()));
				return false;
			}
		}
		else {
			logger.log(Level.FINE, "Has mismatch in repeats");
			ArrayList<Integer> index1 = allIndex(seq.substring(p, q+len), repeat1);
			ArrayList<Integer> index2 = allIndex(seq.substring(p, q+len), repeat2);
			if (index1.size()+index2.size()==2) {
				return true;
			}
			else {
				logger.log(Level.FINE, "More than two repeats:" + Integer.toString(index1.size() + index2.size()));
				return false;
			}
		}
	}

	private int checkReadsForDR (List<PairedEndReads> reads) {
		String seq;
		int withAtLeastTwoDR = 0;
		RepeatInfo repeatinfo = new RepeatInfo();
		for (PairedEndReads r: reads) {
			repeatinfo.reset();
			seq = r.getFirst().getSeq();
			findRepeatsSimple(seq, minRepeatLength, maxRepeatLength, minSpacerLength, 
					maxSpacerLength, repeatinfo);
			if (repeatinfo.isValid() && goodRepeat(repeatinfo.getRepeat())) {
				withAtLeastTwoDR++;
				try {
					this.bufferedWriterDr.write(r.getFirst().getID());
					this.bufferedWriterDr.write('\t');
					this.bufferedWriterDr.write(Integer.toString(repeatinfo.getPositions().get(0)));
					this.bufferedWriterDr.write('/');
					this.bufferedWriterDr.write(Integer.toString(repeatinfo.getPositions().get(1)));
					this.bufferedWriterDr.write('/');
					String repeat1 = repeatinfo.getRepeat();
					String repeat2 = seq.substring(repeatinfo.getPositions().get(1), repeatinfo.getPositions().get(1)+repeat1.length());
					this.bufferedWriterDr.write(repeat1);
					this.bufferedWriterDr.write('/');
					this.bufferedWriterDr.write(repeat2);
					this.bufferedWriterDr.newLine();
					this.bufferedWriterDr.write(repeat1);
					this.bufferedWriterDr.newLine();
					this.bufferedWriterDr.write(repeat2);
					this.bufferedWriterDr.newLine();

					this.bufferedWriterReads.write(r.getFirst().getID());
					this.bufferedWriterReads.newLine();
					this.bufferedWriterReads.write(seq);
					this.bufferedWriterReads.newLine();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			repeatinfo.reset();
			seq = r.getSecond().getSeq();
			findRepeatsSimple(seq, minRepeatLength, maxRepeatLength, minSpacerLength, 
					maxSpacerLength, repeatinfo);
			if (repeatinfo.isValid() && goodRepeat(repeatinfo.getRepeat())) {
				withAtLeastTwoDR++;
				try {
					this.bufferedWriterDr.write(r.getSecond().getID());
					this.bufferedWriterDr.write('\t');
					this.bufferedWriterDr.write(Integer.toString(repeatinfo.getPositions().get(0)));
					this.bufferedWriterDr.write('/');
					this.bufferedWriterDr.write(Integer.toString(repeatinfo.getPositions().get(1)));
					this.bufferedWriterDr.write('/');
					String repeat1 = repeatinfo.getRepeat();
					String repeat2 = seq.substring(repeatinfo.getPositions().get(1), repeatinfo.getPositions().get(1)+repeat1.length());
					this.bufferedWriterDr.write(repeat1);
					this.bufferedWriterDr.write('/');
					this.bufferedWriterDr.write(repeat2);
					this.bufferedWriterDr.newLine();
					this.bufferedWriterDr.write(repeat1);
					this.bufferedWriterDr.newLine();
					this.bufferedWriterDr.write(repeat2);
					this.bufferedWriterDr.newLine();

					this.bufferedWriterReads.write(r.getSecond().getID());
					this.bufferedWriterReads.newLine();
					this.bufferedWriterReads.write(seq);
					this.bufferedWriterReads.newLine();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}
		return withAtLeastTwoDR;
	}

	private void findRepeatsSimple(String sequence, int minRepeatLength, int maxRepeatLength, 
			int minSpacerLength, int maxSpacerLength, RepeatInfo repeatinfo)
	{
		int sequenceLength = sequence.length();
		String pattern;
		String repeat;

		int searchEnd = sequenceLength - maxRepeatLength - maxSpacerLength - minRepeatLength;
		if (searchEnd<=0) {
			searchEnd = sequenceLength - minRepeatLength - minSpacerLength - minRepeatLength;
		}
		for (int j = 0; j <= searchEnd; j = j+1)
		{
			pattern = sequence.substring(j, j + seedSize);
			int patternInTextIndex = sequence.indexOf(pattern, j+minRepeatLength+minSpacerLength);
			if (patternInTextIndex >= 0)
			{   
				int extendLen = extend(sequence, j, patternInTextIndex, seedSize, this.maxMismatch);
				if (seedSize+extendLen < minRepeatLength) {
					continue;
				}
				if (minSpacerLength <= (patternInTextIndex-(j+seedSize+extendLen))) { //spacer length
					if (seedSize+extendLen <= maxRepeatLength) { // repeat length
						repeat = sequence.substring(j, j+seedSize+extendLen);
						Logger.getLogger(RepeatReadsRecuiter.class.getName()).log(Level.FINE, "repeat: " + repeat);
						if (goodSeq(sequence, j, patternInTextIndex, repeat.length())) {
							repeatinfo.setRepeat(repeat);
							repeatinfo.addPosition(j);
							repeatinfo.addPosition(patternInTextIndex);
							repeatinfo.setValid(true);
							return;
						}
					}
					else {
						// if there is a repeat length larger than the max repeat length, we do not need to check others.
						Logger.getLogger(RepeatReadsRecuiter.class.getName()).log(Level.FINE, "Repeat length too large");
						repeatinfo.setValid(false);
						return;
					}
				} else {
					// Found a repeat but the spacer is too small, no need to check others. 
					Logger.getLogger(RepeatReadsRecuiter.class.getName()).log(Level.FINE, "Spacer too small");
					repeatinfo.setValid(false);
					return;
				}
			}
		}
		Logger.getLogger(RepeatReadsRecuiter.class.getName()).log(Level.FINE, "No repeat found");
		repeatinfo.setValid(false);
	}


	private int extend(String seq, int pos1, int pos2, int windowLen, int maxmismatch) {
		int numExtend = 0;
		int numExtendAfter = 0;
		int mismatch = 0;
		int i = pos1 + windowLen;
		int j = pos2 + windowLen;
		//System.out.println(seq.substring(i));
		//System.out.println(seq.substring(j));
		int seqLen = seq.length();
		while (j<seqLen) {
			if (seq.charAt(i) == seq.charAt(j)) {
				i++;
				j++;
				if (mismatch>=1) {
					numExtendAfter++;
				}
				else {
					numExtend++;
				}
			}
			else {
				if (mismatch < maxmismatch) {
					i++;
					j++;
					mismatch++;
					continue;
				}
				else {
					break;
				}
			}
		}
		logger.log(Level.FINE, "Values:");
		logger.log(Level.FINE, Integer.toString(numExtend));
		logger.log(Level.FINE, Integer.toString(mismatch));
		logger.log(Level.FINE, Integer.toString(numExtendAfter));
		if (numExtendAfter <= 3)
			return numExtend;
		else{
			int totalExend = numExtend + mismatch + numExtendAfter;
			if (seq.charAt(pos1+windowLen+totalExend-1) != seq.charAt(pos2+windowLen+totalExend-1)) {
				totalExend = totalExend - 1;
			}
			return totalExend;
		}
	}
}


class ReadFilter implements Runnable{
	private HashSet<String> kmerSet1;
	private HashSet<String> kmerSet2;
	private HashMap<String, String> DRReads;
	private int kmerSize;
	private BufferedWriter bufferedWriter;
	private BufferedWriter bufferedWriterNameMapping;
	private final BlockingQueue <List<PairedEndReads>>  sharedQueue;
	private int internalID;
	private int nFilterThreads;
	private Logger logger = Logger.getLogger(ReadFilter.class.getName());
	public ReadFilter (BlockingQueue <List<PairedEndReads>> sharedQueue, int kmerSize,
			Set<String> kmerSet1, Set<String>kmerSet2, HashMap<String, String> dRReads2, String outName, String outNamemapping, int ID, int nThreads) throws IOException {
		this.sharedQueue = sharedQueue;
		this.kmerSet1 = (HashSet<String>) kmerSet1;
		this.kmerSet2 = (HashSet<String>) kmerSet2;
		this.kmerSize = kmerSize;
		this.bufferedWriter = new BufferedWriter(new FileWriter(new File(outName)));
		this.bufferedWriterNameMapping = new BufferedWriter(new FileWriter(new File(outNamemapping)));
		this.internalID = ID;
		this.nFilterThreads = nThreads;
		this.DRReads = (HashMap<String, String>)dRReads2;
	}

	@Override
	public void run() {
		int readsProcessed = 0;
		int curID = this.internalID;
		StringBuilder sb = new StringBuilder();
		HashMap<String, String> posMap = new HashMap<String, String>();
		while(true){
			try {
				logger.log(Level.INFO, "Taking a new list");
				List<PairedEndReads> inList = sharedQueue.take();
				logger.log(Level.INFO, "New list size: " + Integer.toString(inList.size()));
				readsProcessed += 2 * inList.size();
				if (inList.size() == 0) {
					logger.log(Level.INFO, "Done. ReadFilter thread exit..");
					try {
						bufferedWriter.close();
						bufferedWriterNameMapping.close();
					} catch (IOException e) {
						e.printStackTrace();
					}
					break;
				}
				else {
					posMap.clear();
					List<PairedEndReads> outList = filterReads(inList, posMap);
					float percent = 100 * (float)(outList.size()) / inList.size();
					logger.log(Level.INFO, String.format("Total reads proccessed by this thread:  %d reads, recuited %.3f pecent of the input reads",
							readsProcessed, percent));
					for (PairedEndReads r: outList) {
						try {
							sb.setLength(0);
							sb.append(curID);
							if (this.DRReads.containsKey(r.getFirst().getID())) {
								sb.append(".r.1\t");
								sb.append(this.DRReads.get(r.getFirst().getID()));
							}
							else {
								sb.append(".1");
								if (posMap.containsKey(r.getFirst().getID())) {
									sb.append('\t');
									sb.append(posMap.get(r.getFirst().getID()));
								}
							} 	
							bufferedWriterNameMapping.write(r.getFirst().getID().substring(1));
							bufferedWriterNameMapping.write("\t");
							bufferedWriterNameMapping.write(sb.toString());
							bufferedWriterNameMapping.newLine();

							bufferedWriter.write(">");
							bufferedWriter.write(sb.toString());
							bufferedWriter.newLine();
							bufferedWriter.write(r.getFirst().getSeq());
							bufferedWriter.newLine();

							sb.setLength(0);
							sb.append(curID);
							if (this.DRReads.containsKey(r.getSecond().getID())) {
								sb.append(".r.2\t");
								sb.append(this.DRReads.get(r.getSecond().getID()));
							}
							else {
								sb.append(".2");
								if (posMap.containsKey(r.getSecond().getID())) {
									sb.append('\t');
									sb.append(posMap.get(r.getSecond().getID()));
								}
							}

							bufferedWriterNameMapping.write(r.getSecond().getID().substring(1));
							bufferedWriterNameMapping.write("\t");
							bufferedWriterNameMapping.write(sb.toString());
							bufferedWriterNameMapping.newLine();

							bufferedWriter.write(">");
							bufferedWriter.write(sb.toString());
							bufferedWriter.newLine();
							bufferedWriter.write(r.getSecond().getSeq());
							bufferedWriter.newLine();
							curID = curID + nFilterThreads;
						} catch (IOException e) {
							e.printStackTrace();
						}
					}
				}
			} catch (InterruptedException ex) {
				Logger.getLogger(ReadFilter.class.getName()).log(Level.SEVERE, null, ex);
			}
		}
	}

	private List<PairedEndReads> filterReads (List<PairedEndReads> reads, Map<String, String> posInfo) {
		String seq1, seq2;
		int end;
		boolean firstMatch, secondMatch;
		StringBuilder sb = new StringBuilder();
		List<PairedEndReads> outList= new ArrayList<PairedEndReads>();
		for (PairedEndReads r: reads) {
			sb.setLength(0);
			firstMatch = false;
			seq1 = r.getFirst().getSeq();
			end = seq1.length() - kmerSize;
			for (int i=0; i<=end; i++) {
				String curKmer = seq1.substring(i, i + kmerSize);
				if (kmerSet1.contains(curKmer)) {
					firstMatch = true;
					outList.add(r);
					sb.append(i);
					sb.append("/");
					sb.append(curKmer);
					sb.append("/");
					sb.append(1);
					posInfo.put(r.getFirst().getID(), sb.toString());
					break;
				} 
				else if (kmerSet2.contains(curKmer)) {
					firstMatch = true;
					outList.add(r);
					sb.append(i);
					sb.append("/");
					sb.append(curKmer);
					sb.append("/");
					sb.append(2);
					posInfo.put(r.getFirst().getID(), sb.toString());
					break;
				}
			}
			//if (!firstMatch) {
			sb.setLength(0);
			seq2 = r.getSecond().getSeq();
			end = seq2.length() - kmerSize;
			secondMatch = false;
			for (int i=0; i<=end; i++) {
				String curKmer = seq2.substring(i, i + kmerSize);
				if (kmerSet1.contains(curKmer)) {
					if (!firstMatch ) {
						outList.add(r);
					}
					secondMatch = true;
					sb.append(i);
					sb.append("/");
					sb.append(curKmer);
					sb.append("/");
					sb.append(1);
					posInfo.put(r.getSecond().getID(), sb.toString());
					break;
				}
				else if (kmerSet2.contains(curKmer)) {
					if (! firstMatch) {
						outList.add(r);
					}
					secondMatch = true;
					sb.append(i);
					sb.append("/");
					sb.append(curKmer);
					sb.append("/");
					sb.append(2);
					posInfo.put(r.getSecond().getID(), sb.toString());
					break;
				}
			}
			if (!firstMatch && !secondMatch) {
				String lcs = LongestCommonSubstring.lcs(seq1, ReadRecuiter.reverseCompliment(seq2));
				if (lcs.length() < 20 || lcs.length() > 60)
					continue;
				posInfo.put(r.getFirst().getID(), "p/"+lcs);
				posInfo.put(r.getSecond().getID(), "p/"+lcs);
				outList.add(r);
			}
			//}
		}
		return outList;
	}
}
