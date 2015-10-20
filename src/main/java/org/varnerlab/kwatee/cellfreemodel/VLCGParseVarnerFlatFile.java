package org.varnerlab.kwatee.cellfreemodel;

// import -
import org.varnerlab.kwatee.cellfreemodel.model.VLCGAllostericControlModel;
import org.varnerlab.kwatee.cellfreemodel.model.VLCGMetabolicReactionModel;
import org.varnerlab.kwatee.cellfreemodel.model.VLCGMetaboliteSymbol;
import org.varnerlab.kwatee.cellfreemodel.parserdelegate.VLCGAllostericParserDelegate;
import org.varnerlab.kwatee.cellfreemodel.parserdelegate.VLCGMetabolismParserDelegate;
import org.varnerlab.kwatee.cellfreemodel.parserdelegate.VLCGParserHandlerDelegate;
import org.varnerlab.kwatee.foundation.*;
import org.w3c.dom.Document;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathFactory;
import java.io.*;
import java.util.*;

import org.sbml.libsbml.*;
import org.xml.sax.InputSource;

/**
 * Copyright (c) 2015 Varnerlab,
 * School of Chemical Engineering,
 * Purdue University, West Lafayette IN 46077 USA.
 * <p>
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * <p>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 * <p>
 * Created by jeffreyvarner on 10/8/15.
 */


public class VLCGParseVarnerFlatFile implements VLCGInputHandler {

    // instance variables -
    private VLCGTransformationPropertyTree _transformation_properties_tree = null;
    private XPathFactory _xpath_factory = XPathFactory.newInstance();
    private XPath _xpath_object = _xpath_factory.newXPath();

    private Vector<VLCGMetabolicReactionModel> _reaction_vector = new Vector<VLCGMetabolicReactionModel>();
    private Vector<String> _alphabet_vector = new Vector<String>();
    private Hashtable _resource_table = new Hashtable();
    private Model _model_wrapper = null;

    private Vector<VLCGMetabolicReactionModel> reaction_model_vector = new Vector<VLCGMetabolicReactionModel>();
    private Vector<VLCGAllostericControlModel> allosteric_model_vector = new Vector<VLCGAllostericControlModel>();

    public static final String CELLFREE_METABOLISM_MODEL_TREE = "CELLFREE_METABOLISM_MODEL_TREE";
    public static final String CELLFREE_CONTROL_MODEL_TREE = "CELLFREE_CONTROL_MODEL_TREE";

    public VLCGParseVarnerFlatFile() {
    }


    public void setPropertiesTree(VLCGTransformationPropertyTree properties_tree) {

        if (properties_tree == null){
            return;
        }


        _transformation_properties_tree = properties_tree;
    }

    public Object getResource(Object object) throws Exception {

        // Create a new modelWrapper -
        _model_wrapper = new Model(3,1);

        // Add species to model -
        _addSpeciesToModel();

        // Add the reactions to the model -
        _addReactionsToModel();

        // Create the control function -
        Document control_tree = _createAllostericControlTree();

        // Wrap the control tree -
        VLCGAllostericControlTreeWrapper allosteric_wrapper = new VLCGAllostericControlTreeWrapper();
        allosteric_wrapper.set_allosteric_control_tree(control_tree);

        // package -
        _resource_table.put(CELLFREE_METABOLISM_MODEL_TREE,_model_wrapper);
        _resource_table.put(CELLFREE_CONTROL_MODEL_TREE,allosteric_wrapper);

        // return the model wrapper -
        return _resource_table;
    }

    public void loadResource(Object object) throws Exception {

        // method variables -

        // Where is the file that I need to load?
        String resource_file_path = _transformation_properties_tree.lookupKwateeNetworkFilePath();
        if (resource_file_path != null){

            // ok, we have what appears to be a path, read the VFF file at this location -
            //_readMetabolicReactionData(resource_file_path, reaction_model_vector);
            _readMetabolicModelDescriptionFile(resource_file_path,reaction_model_vector,allosteric_model_vector);

            // We now have the reactions loaded into a vector.
            // We need to check for reversible reactions, and split into two irreversible fragments -
            Iterator<VLCGMetabolicReactionModel> iterator = reaction_model_vector.iterator();
            while (iterator.hasNext()){

                // Grab the reaction model -
                VLCGMetabolicReactionModel reaction_model = iterator.next();

                // is this reaction reversible?
                String strReverseFlag = (String)reaction_model.getMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_REVERSE);
                String strRxnNameTest = (String)reaction_model.getMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_NAME);

                if (strReverseFlag.equalsIgnoreCase("-inf") && strRxnNameTest != "")
                {
                    // Ok if we here, then I have a reversible reaction
                    String strReactants = (String)reaction_model.getMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_REACTANTS);
                    String strProducts = (String)reaction_model.getMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_PRODUCTS);
                    String strRxnName = (String)reaction_model.getMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_NAME);

                    // Create a 2 new records, one with the forward reaction and one for the reverse -
                    VLCGMetabolicReactionModel recForward = new VLCGMetabolicReactionModel();
                    recForward.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_REVERSE, "0");
                    recForward.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_FORWARD,"inf;");
                    recForward.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_REACTANTS,strReactants);
                    recForward.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_PRODUCTS,strProducts);
                    recForward.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_NAME,strRxnName);
                    recForward.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_TYPE,"ON-RATE");
                    recForward.doExecute();

                    VLCGMetabolicReactionModel recReverse = new VLCGMetabolicReactionModel();
                    recReverse.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_REVERSE, "0");
                    recReverse.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_FORWARD, "inf;");
                    recReverse.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_PRODUCTS, strReactants);
                    recReverse.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_REACTANTS, strProducts);
                    recReverse.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_NAME, strRxnName + "_R");
                    recReverse.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_TYPE, "OFF-RATE");
                    recReverse.doExecute();

                    // Store the rxnObject in a tmpVector
                    _reaction_vector.addElement(recForward);
                    _reaction_vector.addElement(recReverse);
                }
                else if (strRxnNameTest!="")
                {
                    // Before I add the record to the rxnObject - give it a type -
                    reaction_model.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_TYPE,"CAT-RATE");
                    reaction_model.doExecute();

                    // Store the rxnObject in a tmpVector
                    _reaction_vector.addElement(reaction_model);
                }
            } // end while

            // create a new alphabet vector -
            _createUniqueAlphabet();

            // Sort according to order file
            _orderMyUniqueMetaboliteAlphabet();
        }
        else {
            throw new Exception("ERROR: Missing resource file path. Can't find VFF to parse.");
        }
    }


    private Document _createAllostericControlTree() throws Exception {

        // Method variables -
        StringBuffer xml_buffer = new StringBuffer();
        DocumentBuilder document_builder = null;
        Document control_tree = null;

        // Populate the string buffer -
        xml_buffer.append("<?xml version=\"1.0\" encoding=\"utf-8\" standalone=\"yes\"?>\n");
        xml_buffer.append("<ControlModel>\n");
        xml_buffer.append("\t<allosteric>\n");

        // process the list of allosteric connections -
        int counter = 1;
        Iterator<VLCGAllostericControlModel> allosteric_iterator = allosteric_model_vector.iterator();
        while (allosteric_iterator.hasNext()){

            // Get the control model -
            VLCGAllostericControlModel control_model = allosteric_iterator.next();

            // Get the data from the model -
            String actor = (String)control_model.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_ACTOR);
            String type = (String)control_model.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_TYPE);
            String target = (String)control_model.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_TARGET);
            String name = (String)control_model.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_NAME);

            // write the record -
            xml_buffer.append("\t\t<connection name=\"");
            xml_buffer.append(name);
            xml_buffer.append("\" type=\"");
            xml_buffer.append(type);
            xml_buffer.append("\" actor=\"");
            xml_buffer.append(actor);
            xml_buffer.append("\" target=\"");
            xml_buffer.append(target);
            xml_buffer.append("\"/>\n");
        }

        xml_buffer.append("\t</allosteric>\n");
        xml_buffer.append("</ControlModel>\n");

        System.out.println(xml_buffer.toString());

        // Convert the string buffer into an XML Document object -
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        document_builder = factory.newDocumentBuilder();
        control_tree = document_builder.parse(new InputSource(new StringReader(xml_buffer.toString())));

        // return the control tree -
        return control_tree;
    }

    private void _orderMyUniqueMetaboliteAlphabet() throws Exception {

        // method variables -
        Vector<String> species_order_vector = new Vector<String>();

        // Get the path to the order file -
        String path_to_order_file = _transformation_properties_tree.lookupKwateeSpeciesOrderFilePath();
        if (path_to_order_file == null){
            throw new Exception("Path to the species order file is null?");
        }

        // Is there a file at the end of this rainbow?
        File order_file = new File(path_to_order_file);
        if (order_file.exists() && !order_file.isDirectory()) {

            // ok - we have the order file, load that data into a temp vector -
            BufferedReader inReader = new BufferedReader(new FileReader(order_file));
            String record = null;
            while ((record = inReader.readLine()) != null) {

                species_order_vector.addElement(record);
            }

            // close the reader -
            inReader.close();

            // ok, we have the species loaded, let's reorder my alphabet vector -
            if (species_order_vector.size()>0){

                // iterate through -
                ListIterator order_list_iterator = species_order_vector.listIterator();
                while (order_list_iterator.hasNext()){

                    // Get the species -
                    String species_symbol = (String)order_list_iterator.next();

                    // Is this symbol in the alphabet?
                    if (_alphabet_vector.contains(species_symbol)){

                        // ok, remove the symbol from the alphabet vector, and put at the end -
                        _alphabet_vector.remove(species_symbol);
                        _alphabet_vector.addElement(species_symbol);
                    }
                }
            }
        }
        else {
            throw new Exception("File at path "+path_to_order_file+" was not found?");
        }
    }

    private void _createUniqueAlphabet() throws Exception {

        // Get the iterator of reaction objects
        Iterator<VLCGMetabolicReactionModel> iterReactions = _reaction_vector.iterator();
        while (iterReactions.hasNext())
        {
            // Ok, get the reactants and products -
            VLCGMetabolicReactionModel rxnObject = iterReactions.next();
            Iterator<VLCGMetaboliteSymbol> iterReactants = rxnObject.getReactants();
            while (iterReactants.hasNext())
            {
                // Get StateSymbol object
                VLCGMetaboliteSymbol symbolReactants = iterReactants.next();

                String testSymbol = symbolReactants.getSymbol();
                if (!_alphabet_vector.contains(testSymbol) && !testSymbol.equalsIgnoreCase("[]"))
                {
                    _alphabet_vector.addElement(testSymbol);
                }
            }

            Iterator<VLCGMetaboliteSymbol> iterProducts = rxnObject.getProducts();
            while (iterProducts.hasNext())
            {
                // Get StateSymbol object
                VLCGMetaboliteSymbol symbolProducts = iterProducts.next();

                String testSymbol = symbolProducts.getSymbol();
                if (!_alphabet_vector.contains(testSymbol) && !testSymbol.equalsIgnoreCase("[]"))
                {
                    _alphabet_vector.addElement(testSymbol);
                }
            }
        }
    }

    private void _readMetabolicModelDescriptionFile(String fileName,
                                                    Vector<VLCGMetabolicReactionModel> reaction_vector,
                                                    Vector<VLCGAllostericControlModel> allosteric_vector) throws Exception {

        // method allocation -
        VLCGParserHandlerDelegate parser_delegate = null;

        // check -
        if (fileName == null ||
                reaction_vector == null ||
                allosteric_vector == null){

            return;
        }

        BufferedReader inReader = new BufferedReader(new FileReader(fileName));
        String dataRecord=null;
        while ((dataRecord=inReader.readLine()) != null) {

            int whitespace = dataRecord.length();

            // Need to check to make sure I have do not have a comment
            if (!dataRecord.contains("//") && whitespace != 0) {

                // is this a handler directive?
                if (dataRecord.equalsIgnoreCase("[kwatee.metabolic.reaction.handler]") == true){

                    // fire up the reaction handler for metabolism -
                    parser_delegate = new VLCGMetabolismParserDelegate();
                }
                else if(dataRecord.equalsIgnoreCase("[kwatee.allosteric.control.handler]") == true){

                    // fire up the allosteric handler -
                    parser_delegate = new VLCGAllostericParserDelegate();
                }

                if (!dataRecord.contains("[kwatee.metabolic.reaction.handler]") &&
                        !dataRecord.equalsIgnoreCase("[kwatee.allosteric.control.handler]")) {

                    //System.out.println(" WTF - "+dataRecord+" test = "+(!dataRecord.contains("[kwatee.metabolic.reaction.handler]")));

                    // build the parsed object -
                    Object parsed_object = parser_delegate.parseLine(dataRecord);

                    // Check - is the parsed_object a reaction model, or a control model -
                    if (parsed_object instanceof VLCGMetabolicReactionModel){
                        reaction_vector.addElement((VLCGMetabolicReactionModel) parsed_object);
                        //System.out.println("Adding ...");
                    }
                    else {
                        allosteric_vector.addElement((VLCGAllostericControlModel) parsed_object);
                    }
                }
            }
        }

        // close -
        inReader.close();
    }

    private void _readMetabolicReactionData(String fileName,Vector<VLCGMetabolicReactionModel> vector) throws Exception {
        // check -
        if (fileName == null || vector == null){
            return;
        }

        // Ok, now down to bizness...
        BufferedReader inReader = new BufferedReader(new FileReader(fileName));
        String dataRecord=null;

        while ((dataRecord=inReader.readLine())!=null)
        {
            // When I get here I have a data record, I need to
            // parse the record
            int whitespace = dataRecord.length();

            // Need to check to make sure I have do not have a comment
            if (!dataRecord.contains("//") && whitespace !=0) {

                // Create a data record wrapper
                VLCGMetabolicReactionModel metabolic_reaction_model = new VLCGMetabolicReactionModel();

                StringTokenizer tokenizer=new StringTokenizer(cleanRecord(dataRecord),",",false);
                int intCounter = 1;
                while (tokenizer.hasMoreElements())
                {

                    // Get a data from the tokenizer -
                    Object dataChunk = tokenizer.nextToken();

                    if (intCounter==1) {
                        metabolic_reaction_model.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_NAME, dataChunk);
                    }
                    else if (intCounter==2) {
                        String strTmp = ((String)dataChunk).replace("-", "_");
                        metabolic_reaction_model.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_REACTANTS, strTmp);
                    }
                    else if (intCounter==3) {
                        String strTmp = ((String)dataChunk).replace("-", "_");
                        metabolic_reaction_model.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_PRODUCTS, strTmp);
                    }
                    else if (intCounter==4) {
                        metabolic_reaction_model.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_REVERSE, dataChunk);
                    }
                    else if (intCounter==5) {
                        metabolic_reaction_model.setMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_REVERSE, dataChunk);
                    }

                    // increment the counter
                    intCounter++;
                }

                // Add the record to the vector
                vector.addElement(metabolic_reaction_model);
            }
        }
    }

    private String cleanRecord(String record) {
        String rString = record.replaceAll(",,",",NULL,");
        return(rString);
    }

    private void _addReactionsToModel() throws Exception {
        // Method attributes -
        int counter = 0;

        // Get the iterator of FFReactionObjects -
        Iterator reaction_iterator = _reaction_vector.iterator();

        // Iterate through the list of VFF reactions and populate the SBML
        while (reaction_iterator.hasNext()) {
            // Get the ffRxnObj -
            VLCGMetabolicReactionModel ffRxnObj = (VLCGMetabolicReactionModel)reaction_iterator.next();

            // Create a new SBML reaction object -
            Reaction sbml_reaction = _model_wrapper.createReaction();

            // Get the reaction name?
            // String reaction_name = (String)ffRxnObj.getMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_NAME);

            // flat files always are irreversible -
            sbml_reaction.setReversible(false);
            // sbml_reaction.setName(reaction_name);

            // Process the reactants -
            Iterator iterator_reactants = ffRxnObj.getReactants();
            while (iterator_reactants.hasNext()) {

                // Get the state symbol -
                VLCGMetaboliteSymbol symbol = (VLCGMetaboliteSymbol)iterator_reactants.next();

                // Get the symbol -
                String strSymbol = symbol.getSymbol();

                // Create and configure species reference -
                SpeciesReference specRef = sbml_reaction.createReactant();

                // Check to see if we have a [] -
                if (strSymbol.equalsIgnoreCase("[]")) {
                    specRef.setSpecies("[]");
                }
                else {
                    specRef.setSpecies(strSymbol);
                }

                // We need to check VFF puts negatives in -
                double tmpCoeff = Math.abs(ffRxnObj.getCoefficientValue(strSymbol));
                specRef.setStoichiometry(tmpCoeff);

                // Add the configured species reference -
                sbml_reaction.addReactant(specRef);
            }

            // Configure the reaction rate - add the products -
            Iterator iter_products = ffRxnObj.getProducts();
            while (iter_products.hasNext()) {
                // Get the state symbol -
                VLCGMetaboliteSymbol symbol = (VLCGMetaboliteSymbol)iter_products.next();

                // Create a new SpeciesRef -
                SpeciesReference specRef = sbml_reaction.createProduct();
                specRef.setSpecies(symbol.getSymbol());
                specRef.setStoichiometry(ffRxnObj.getCoefficientValue(symbol.getSymbol()));

                // Add the species ref -
                sbml_reaction.addProduct(specRef);
            }

            // Get the reaction string and add as the name -
            sbml_reaction.setName(ffRxnObj.getReactionString());
            sbml_reaction.setId((String)ffRxnObj.getMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_NAME));

            // add the reaction to the model -
            _model_wrapper.addReaction(sbml_reaction);

            // update the counter -
            counter++;
        }

        // Add enzyme degrdation rates to model -
        reaction_iterator = _reaction_vector.iterator();
        int reaction_counter = 1;
        while (reaction_iterator.hasNext()){

            // Each reaction has a enzyme which degrades -
            VLCGMetabolicReactionModel reaction_model = (VLCGMetabolicReactionModel)reaction_iterator.next();

            // Create a new SBML reaction object -
            Reaction sbml_reaction = _model_wrapper.createReaction();
            sbml_reaction.setReversible(false);

            // Create and configure species reference -
            SpeciesReference specRef = sbml_reaction.createReactant();
            String enzyme = "E_"+reaction_model.getMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_NAME);
            specRef.setSpecies(enzyme);
            specRef.setStoichiometry(1.0);

            // Add the configured species reference -
            sbml_reaction.addReactant(specRef);

            // Add reaction to model -
            sbml_reaction.setName("Degradation: "+enzyme+" = []");
            sbml_reaction.setId(enzyme);
            _model_wrapper.addReaction(sbml_reaction);

            // update counter -
            reaction_counter++;
        }
    }

    private void _addSpeciesToModel() throws Exception {

        // Method attributes -
        int counter = 0;

        // Ok, get the iterator for species so I can add these to the model -
        Iterator species_iter = _alphabet_vector.iterator();

        // Process each species -
        while (species_iter.hasNext()){

            // Get the species -
            String spName = (String)species_iter.next();

            // We need to check for -
            spName = spName.replace('-', '_');

            // Configure the species object -
            Species newSpecies = _model_wrapper.createSpecies();
            newSpecies.setName(spName);
            newSpecies.setId(spName);

            // Get the model name -
            String strModelName = _transformation_properties_tree.lookupKwateeModelName();
            newSpecies.setCompartment(strModelName);

            // if the initial condition is there, use it
            // we are assuming the ic file is in same order as this
            newSpecies.setInitialConcentration(0.0);

            // Add the configured species to the model -
            _model_wrapper.addSpecies(newSpecies);
            counter++;
        }

        // for a CF model, we also have enzymes for each of the reactions -
        Iterator reaction_iterator = _reaction_vector.iterator();
        int reaction_counter = 1;
        while (reaction_iterator.hasNext()){

            // Reaction object -
            VLCGMetabolicReactionModel reaction = (VLCGMetabolicReactionModel)reaction_iterator.next();

            // Name and id -
            String name = "E_"+reaction.getMetabolicReactionComponent(VLCGMetabolicReactionModel.REACTION_NAME);

            // Make a species -
            Species newSpecies = _model_wrapper.createSpecies();
            newSpecies.setName(name);
            newSpecies.setId(name);

            // Get the model name -
            String strModelName = _transformation_properties_tree.lookupKwateeModelName();
            newSpecies.setCompartment(strModelName);

            // CF enzymes always start at 1 -
            newSpecies.setInitialConcentration(1.0);

            // Add the configured species to the model -
            _model_wrapper.addSpecies(newSpecies);
            reaction_counter++;
        }
    }


    /**
     * Loads the SWIG-generated libSBML Java module when this class is
     * loaded, or reports a sensible diagnostic message about why it failed.
     */
    static
    {
        try
        {
            System.loadLibrary("sbmlj");
            // For extra safety, check that the jar file is in the classpath.
            Class.forName("org.sbml.libsbml.libsbml");
        }
        catch (UnsatisfiedLinkError e)
        {
            System.err.println("Error encountered while attempting to load libSBML:");
            System.err.println("Please check the value of your "
                    + (System.getProperty("os.name").startsWith("Mac OS")
                    ? "DYLD_LIBRARY_PATH" : "LD_LIBRARY_PATH") +
                    " environment variable and/or your" +
                    " 'java.library.path' system property (depending on" +
                    " which one you are using) to make sure it list the" +
                    " directories needed to find the " +
                    System.mapLibraryName("sbmlj") + " library file and" +
                    " libraries it depends upon (e.g., the XML parser).");
            System.exit(1);
        }
        catch (ClassNotFoundException e)
        {
            System.err.println("Error: unable to load the file 'libsbmlj.jar'." +
                    " It is likely that your -classpath command line " +
                    " setting or your CLASSPATH environment variable " +
                    " do not include the file 'libsbmlj.jar'.");
            e.printStackTrace();

            System.exit(1);
        }
        catch (SecurityException e)
        {
            System.err.println("Error encountered while attempting to load libSBML:");
            e.printStackTrace();
            System.err.println("Could not load the libSBML library files due to a"+
                    " security exception.\n");
            System.exit(1);
        }
    }
}
