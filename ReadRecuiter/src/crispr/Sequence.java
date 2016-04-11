package crispr;


public class Sequence {
	public Sequence(String id, String seq) {
		this.identifier = id;
		this.seq = seq;
	}

	public String getID() {
		return identifier;
	}

	public String getSeq() {
		return seq;
	}

	private String identifier;
	private String seq;
}
