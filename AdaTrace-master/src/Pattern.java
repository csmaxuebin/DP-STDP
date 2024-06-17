import java.util.List;


public class Pattern {
	
	private List<Cell> cellsInPattern;
	
	public Pattern (List<Cell> elts) {
		this.cellsInPattern = elts;
	}
	
	public boolean isEqual(Pattern p2) {
		if (this.cellsInPattern.size() != p2.cellsInPattern.size())
			return false;
		for (int i = 0; i < this.cellsInPattern.size(); i++) {
			if (this.cellsInPattern.get(i) != p2.cellsInPattern.get(i))
				return false;
		}
		return true;
	}
	
	@Override 
	public boolean equals (Object other){
	    if (other == null) return false;
	    if (other == this) return true;
	    if (!(other instanceof Pattern))return false;//如果 object 是 class 的一个实例，则 instanceof 运算符返回 true。如果 object 不是指定类的一个实例，或者 object 是 null，则返回 false
	    Pattern otherMyClass = (Pattern)other;
	    if (this.isEqual(otherMyClass))
	    	return true;
	    else
	    	return false;
	}
	
	@Override
    public int hashCode() {
		final int prime = 31;
		int result = 1;
		for( Cell s : cellsInPattern ){
		    result = result * prime + s.hashCode();
		}
		return result;
    }
	
	@Override
	public String toString() {
		String tbr = "[";
		for (Cell c: cellsInPattern) {
			tbr = tbr + "(" + c.getName() + "),";
		}
		tbr = tbr.substring(0, tbr.length()-1);
		tbr = tbr + "]";
		return tbr;
	}

}
