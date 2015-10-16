package org.varnerlab.kwatee.cellfreemodel;

// imports -
import org.varnerlab.kwatee.foundation.VLCGOutputHandler;
import org.varnerlab.kwatee.foundation.VLCGTransformationPropertyTree;
import org.sbml.libsbml.*;
import org.w3c.dom.Document;

import javax.print.Doc;
import java.util.Hashtable;
import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;

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
 * Created by jeffreyvarner on 10/9/15.
 */

public class VLCGWriteJuliaCellFreeModel implements VLCGOutputHandler {

    // instance variables -
    private VLCGTransformationPropertyTree _transformation_properties_tree = null;
    private VLCGJuliaCellFreeModelDelegate _delegate_object = new VLCGJuliaCellFreeModelDelegate();

    public VLCGWriteJuliaCellFreeModel() {
    }

    /**
     * Generates and writes OutputHandler specific code for a given model.
     * The given model is passed in as an SBML model object.
     * @param object an SBML model object that represents the current model
     * @throws Exception
     */
    public void writeResource(Object object) throws Exception {

        // ok, object should be SMBL -
        Hashtable resource_table = (Hashtable)object;

        // Get SBML model tree -
        Model model_wrapper = (Model)resource_table.get(VLCGParseVarnerFlatFile.CELLFREE_METABOLISM_MODEL_TREE);
        VLCGAllostericControlTreeWrapper control_tree = (VLCGAllostericControlTreeWrapper)resource_table.get(VLCGParseVarnerFlatFile.CELLFREE_CONTROL_MODEL_TREE);

        // Build the data dictionary -
        String fully_qualified_data_dictionary_path = _transformation_properties_tree.lookupKwateeDataDictionaryFilePath();
        String data_dictionary = _delegate_object.buildDataDictionaryBuffer(model_wrapper,control_tree,_transformation_properties_tree);
        write(fully_qualified_data_dictionary_path,data_dictionary);

        // Build the driver -
        String fully_qualified_driver_path = _transformation_properties_tree.lookupKwateeDriverFunctionFilePath();
        String driver_buffer = _delegate_object.buildDriverFunctionBuffer(model_wrapper,_transformation_properties_tree);
        write(fully_qualified_driver_path,driver_buffer);

        // Ok, lets build the stoichiometric matrix -
        int NUMBER_OF_SPECIES = (int)model_wrapper.getNumSpecies();
        int NUMBER_OF_RATES = (int)model_wrapper.getNumReactions();

        // Initialize the stoichiometric matrix -
        String fully_qualified_stoichiometric_matrix_path = _transformation_properties_tree.lookupKwateeStoichiometricMatrixFilePath();
        double[][] stoichiometric_matrix = new double[NUMBER_OF_SPECIES][NUMBER_OF_RATES];
        _delegate_object.buildStoichiometricMatrix(stoichiometric_matrix,model_wrapper);
        writeStoichiometricMatrixToDisk(fully_qualified_stoichiometric_matrix_path,stoichiometric_matrix,model_wrapper);

        // Build the balance equations -
        String fully_qualified_balance_path = _transformation_properties_tree.lookupKwateeBalanceFunctionFilePath();
        String balance_buffer = _delegate_object.buildBalanceFunctionBuffer(model_wrapper,_transformation_properties_tree);
        write(fully_qualified_balance_path,balance_buffer);

        // Build the kinetics equations -
        String fully_qualified_kinetics_path = _transformation_properties_tree.lookupKwateeKineticsFunctionFilePath();
        String kinetics_buffer = _delegate_object.buildKineticsFunctionBuffer(model_wrapper,_transformation_properties_tree);
        write(fully_qualified_kinetics_path,kinetics_buffer);

        // Build the allosteric control equations -
        String fully_qualified_allosteric_control_path = _transformation_properties_tree.lookupKwateeControlFunctionFilePath();
        String control_buffer = _delegate_object.buildAllostericControlFunctionBuffer(model_wrapper,control_tree,_transformation_properties_tree);
        write(fully_qualified_allosteric_control_path,control_buffer);
    }



    /**
     * Used by Kwatee to pass required properties to the InputHandler.
     * When constructing a model there are many properties
     * that are required; these properties may change from one instance to another.
     * For example, the input file path is obviously required but is also variable.
     * Universal is written to store and communicate these properties in an
     * XMLPropTree object.  The XMLPropTree object is then made accessible
     * to the InputHandler via the setProperties Method.
     *
     * When creating a custom plugin, write the setProperties method
     * to properly make use of the "instance specific" information in the
     * XMLPropTree object.
     *
     *
     * For more see the javadoc associated with XMLPropTree class.
     * @param properties_tree an VLCGTransformationPropertyTree object containing various properties
     */
    public void setPropertiesTree(VLCGTransformationPropertyTree properties_tree) {

        if (properties_tree == null){
            return;
        }

        _transformation_properties_tree = properties_tree;
    }

    public Object getResource(Object object) throws Exception {
        return null;
    }

    // Private -
    private void writeStoichiometricMatrixToDisk(String path,double[][] stoichiometric_matrix,Model model_wrapper) throws Exception {

        // Method attribute -
        StringBuffer buffer = new StringBuffer();

        // Get the system dimension -
        int NUMBER_OF_SPECIES = (int)model_wrapper.getNumSpecies();
        int NUMBER_OF_RATES = (int)model_wrapper.getNumReactions();


        for (int scounter=0;scounter<NUMBER_OF_SPECIES;scounter++) {
            for (int rcounter=0;rcounter<NUMBER_OF_RATES;rcounter++) {
                buffer.append(stoichiometric_matrix[scounter][rcounter]);
                buffer.append("\t");
            }
            buffer.append("\n");
        }

        File oFile = new File(path);
        BufferedWriter writer = new BufferedWriter(new FileWriter(oFile));

        // Write buffer to file system and close writer
        writer.write(buffer.toString());
        writer.close();
    }

    private void write(String path,String buffer) throws Exception {

        // Create writer
        File oFile = new File(path);
        BufferedWriter writer = new BufferedWriter(new FileWriter(oFile));

        // Write buffer to file system and close writer
        writer.write(buffer);
        writer.close();
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
