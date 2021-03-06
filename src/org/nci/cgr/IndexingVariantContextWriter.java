package org.nci.cgr;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.LocationAware;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.tribble.index.DynamicIndexCreator;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexCreator;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.samtools.util.PositionalOutputStream;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.nio.file.Path;
public abstract class IndexingVariantContextWriter {
	 private final String name;
	    private final Path location;
	    private final SAMSequenceDictionary refDict;

	    private OutputStream outputStream;
	    private LocationAware locationSource = null;
	    private IndexCreator indexer = null;

	    private IndexingVariantContextWriter(final String name, final Path location, final OutputStream output, final SAMSequenceDictionary refDict) {
	        this.name = name;
	        this.location = location;
	        this.outputStream = output;
	        this.refDict = refDict;
	    }

	    static String DEFAULT_READER_NAME = "Reader Name";

	    /**
	     * Create a VariantContextWriter with an associated index using the default index creator
	     *
	     * @param name  the name of this writer (i.e. the file name or stream)
	     * @param location  the path to the output file
	     * @param output    the output stream to write to
	     * @param refDict   the reference dictionary
	     * @param enableOnTheFlyIndexing    is OTF indexing enabled?
	     */
	    protected IndexingVariantContextWriter(final String name, final File location, final OutputStream output, final SAMSequenceDictionary refDict,
	                                           final boolean enableOnTheFlyIndexing) {
	        this(name, IOUtil.toPath(location), output, refDict, enableOnTheFlyIndexing);
	    }

	    /**
	     * Create a VariantContextWriter with an associated index using the default index creator
	     *
	     * @param name  the name of this writer (i.e. the file name or stream)
	     * @param location  the path to the output file
	     * @param output    the output stream to write to
	     * @param refDict   the reference dictionary
	     * @param enableOnTheFlyIndexing    is OTF indexing enabled?
	     */
	    protected IndexingVariantContextWriter(final String name, final Path location, final OutputStream output, final SAMSequenceDictionary refDict,
	        final boolean enableOnTheFlyIndexing) {
	        this(name, location, output, refDict);

	        if ( enableOnTheFlyIndexing ) {
	            initIndexingWriter(new DynamicIndexCreator(location, IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME));
	        }
	    }

	    /**
	     * Create a VariantContextWriter with an associated index using a custom index creator
	     *
	     * @param name  the name of this writer (i.e. the file name or stream)
	     * @param location  the path to the output file
	     * @param output    the output stream to write to
	     * @param refDict   the reference dictionary
	     * @param enableOnTheFlyIndexing    is OTF indexing enabled?
	     * @param idxCreator    the custom index creator.  NOTE: must be initialized
	     */
	    protected IndexingVariantContextWriter(final String name, final File location, final OutputStream output, final SAMSequenceDictionary refDict,
	                                           final boolean enableOnTheFlyIndexing, final IndexCreator idxCreator) {
	        this(name, IOUtil.toPath(location), output, refDict, enableOnTheFlyIndexing, idxCreator);
	    }

	    /**
	     * Create a VariantContextWriter with an associated index using a custom index creator
	     *
	     * @param name  the name of this writer (i.e. the file name or stream)
	     * @param location  the path to the output file
	     * @param output    the output stream to write to
	     * @param refDict   the reference dictionary
	     * @param enableOnTheFlyIndexing    is OTF indexing enabled?
	     * @param idxCreator    the custom index creator.  NOTE: must be initialized
	     */
	    protected IndexingVariantContextWriter(final String name, final Path location, final OutputStream output, final SAMSequenceDictionary refDict,
	        final boolean enableOnTheFlyIndexing, final IndexCreator idxCreator) {
	        this(name, location, output, refDict);

	        if ( enableOnTheFlyIndexing ) {
	            // TODO: Handle non-Tribble IndexCreators
	            initIndexingWriter(idxCreator);
	        }
	    }

	    private void initIndexingWriter(final IndexCreator idxCreator) {
	        indexer = idxCreator;
	        if (outputStream instanceof LocationAware) {
	            locationSource = (LocationAware)outputStream;
	        } else {
	            final PositionalOutputStream positionalOutputStream = new PositionalOutputStream(outputStream);
	            locationSource = positionalOutputStream;
	            outputStream = positionalOutputStream;
	        }
	    }
	    
	    /** return true is the underlying stream is a PrintStream and 
	     * its checkError returned true. Used to stop linux pipelines 
	     */
	    public boolean checkError() {
	        return (getOutputStream() instanceof PrintStream) && 
	                PrintStream.class.cast(getOutputStream()).checkError();
	    }
	    
	    public OutputStream getOutputStream() {
	        return outputStream;
	    }

	    public String getStreamName() {
	        return name;
	    }

	    public abstract void writeHeader(VCFHeader header);

	    /**
	     * attempt to close the VCF file
	     */
	    public void close() {
	        try {
	            // close the underlying output stream
	            outputStream.close();

	            // close the index stream (keep it separate to help debugging efforts)
	            if (indexer != null) {
	                indexer.setIndexSequenceDictionary(refDict);
	                final Index index = indexer.finalizeIndex(locationSource.getPosition());
	                index.writeBasedOnFeaturePath(location);
	            }


	        } catch (final IOException e) {
	            throw new RuntimeIOException("Unable to close index for " + getStreamName(), e);
	        }
	    }

	    /**
	     * @return the reference sequence dictionary used for the variant contexts being written
	     */
	    public SAMSequenceDictionary getRefDict() {
	        return refDict;
	    }

	    /**
	     * add a record to the file
	     *
	     * @param vc      the Variant Context object
	     */
	    public void add(final VariantContext vc) {
	        // if we are doing on the fly indexing, add the record ***before*** we write any bytes
	        if ( indexer != null )
	            indexer.addFeature(vc, locationSource.getPosition());
	    }

	    /**
	     * Returns a reasonable "name" for this writer, to display to the user if something goes wrong
	     *
	     * @param location
	     * @param stream
	     * @return
	     */
	    protected static final String writerName(final File location, final OutputStream stream) {
	        return writerName(IOUtil.toPath(location), stream);
	    }

	    /**
	     * Returns a reasonable "name" for this writer, to display to the user if something goes wrong
	     *
	     * @param location
	     * @param stream
	     * @return
	     */
	    protected static final String writerName(final Path location, final OutputStream stream) {
	        return location == null ? stream == null ? DEFAULT_READER_NAME : stream.toString() : location.toAbsolutePath().toUri().toString();
	}
}
