package Entity;

/**
 * Entity class for storing gene information retrieved from external APIs
 */
public class GeneInfo {
    private String symbol;
    private String ensemblId;
    private String entrezId;
    private String description;
    private String chromosome;
    private String location;
    private String biotype;
    private String strand;
    private long start;
    private long end;
    private String species;

    public GeneInfo() {
    }

    public GeneInfo(String symbol) {
        this.symbol = symbol;
    }

    // Getters and Setters
    public String getSymbol() {
        return symbol;
    }

    public void setSymbol(String symbol) {
        this.symbol = symbol;
    }

    public String getEnsemblId() {
        return ensemblId;
    }

    public void setEnsemblId(String ensemblId) {
        this.ensemblId = ensemblId;
    }

    public String getEntrezId() {
        return entrezId;
    }

    public void setEntrezId(String entrezId) {
        this.entrezId = entrezId;
    }

    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getChromosome() {
        return chromosome;
    }

    public void setChromosome(String chromosome) {
        this.chromosome = chromosome;
    }

    public String getLocation() {
        return location;
    }

    public void setLocation(String location) {
        this.location = location;
    }

    public String getBiotype() {
        return biotype;
    }

    public void setBiotype(String biotype) {
        this.biotype = biotype;
    }

    public String getStrand() {
        return strand;
    }

    public void setStrand(String strand) {
        this.strand = strand;
    }

    public long getStart() {
        return start;
    }

    public void setStart(long start) {
        this.start = start;
    }

    public long getEnd() {
        return end;
    }

    public void setEnd(long end) {
        this.end = end;
    }

    public String getSpecies() {
        return species;
    }

    public void setSpecies(String species) {
        this.species = species;
    }

    @Override
    public String toString() {
        return "GeneInfo{" +
                "symbol='" + symbol + '\'' +
                ", ensemblId='" + ensemblId + '\'' +
                ", chromosome='" + chromosome + '\'' +
                ", location='" + location + '\'' +
                '}';
    }
}
