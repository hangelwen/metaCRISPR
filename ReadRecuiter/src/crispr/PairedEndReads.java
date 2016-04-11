package crispr;

public class PairedEndReads {
	private Sequence read1;
	private Sequence read2;

	public PairedEndReads(Sequence first, Sequence second) {
		this.read1 = first;
		this.read2 = second;
	}

	public PairedEndReads(String id1, String seq1, String id2, String seq2) {
		this.read1 = new Sequence(id1, seq1);
		this.read2 = new Sequence(id2, seq2);
	}

	public Sequence getFirst() {
		return read1;
	}

	public Sequence getSecond() {
		return read2;
	}
}