package org.varnerlab.kwatee.cellfreemodel;

// import -
import org.sbml.libsbml.*;
import org.varnerlab.kwatee.cellfreemodel.model.VLCGAllostericControlModel;
import org.varnerlab.kwatee.foundation.VLCGCopyrightFactory;
import org.varnerlab.kwatee.foundation.VLCGTransformationPropertyTree;

import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Iterator;
import java.util.Vector;

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


public class VLCGJuliaCellFreeModelDelegate {

    // instance variables -
    private VLCGCopyrightFactory copyrightFactory = VLCGCopyrightFactory.getSharedInstance();
    private java.util.Date today = Calendar.getInstance().getTime();
    private SimpleDateFormat date_formatter = new SimpleDateFormat("MM-dd-yyyy HH:mm:ss");

    public VLCGJuliaCellFreeModelDelegate() {

        // initialize -
        init();
    }

    private void init(){
    }

    public String buildKineticsFunctionBuffer(Model model_tree, VLCGTransformationPropertyTree property_tree) throws Exception {

        // Method variables -
        StringBuffer buffer = new StringBuffer();

        // Copyright notice -
        String copyright = copyrightFactory.getJuliaCopyrightHeader();
        buffer.append(copyright);

        // Get the kinetics function name -
        String kinetics_function_name = property_tree.lookupKwateeKineticsFunctionName();

        // Propulate the buffer -
        buffer.append("function ");
        buffer.append(kinetics_function_name);
        buffer.append("(t,x,data_dictionary)\n");
        buffer.append("# --------------------------------------------------------------------- #\n");
        buffer.append("# ");
        buffer.append(kinetics_function_name);
        buffer.append(".jl was generated using the Kwatee code generation system.\n");
        buffer.append("# Username: ");
        buffer.append(property_tree.lookupKwateeModelUsername());
        buffer.append("\n");
        buffer.append("# Type: ");
        buffer.append(property_tree.lookupKwateeModelType());
        buffer.append("\n");
        buffer.append("# Version: ");
        buffer.append(property_tree.lookupKwateeModelVersion());
        buffer.append("\n");
        buffer.append("# Generation timestamp: ");
        buffer.append(date_formatter.format(today));
        buffer.append("\n");
        buffer.append("# \n");
        buffer.append("# Input arguments: \n");
        buffer.append("# t  - current time \n");
        buffer.append("# x  - state vector \n");
        buffer.append("# data_dictionary - parameter vector \n");
        buffer.append("# \n");
        buffer.append("# Return arguments: \n");
        buffer.append("# rate_vector - rate vector \n");
        buffer.append("# --------------------------------------------------------------------- #\n");
        buffer.append("# \n");
        buffer.append("# Alias the species vector - \n");
        ListOfSpecies listOfSpecies = model_tree.getListOfSpecies();
        long number_of_species = listOfSpecies.size();
        Vector<String> tmp_species_symbol_vector = new Vector<String>();
        for (long species_index = 0;species_index<number_of_species;species_index++){

            // Get the symbol -
            Species species_object = listOfSpecies.get(species_index);
            String species_symbol = species_object.getId();

            // write the symbol =
            //buffer.append("const ");
            buffer.append(species_symbol);
            buffer.append(" = x[");
            buffer.append(species_index + 1);
            buffer.append("];\n");

            // store the symbol for later -
            tmp_species_symbol_vector.add(species_symbol);
        }

        buffer.append("\n");
        buffer.append("# Formulate the kinetic rate vector - \n");
        buffer.append("rate_constant_array = data_dictionary[\"RATE_CONSTANT_ARRAY\"];\n");
        buffer.append("saturation_constant_array = data_dictionary[\"SATURATION_CONSTANT_ARRAY\"];\n");
        buffer.append("rate_vector = Float64[];\n");
        buffer.append("\n");
        ListOfReactions listOfReactions = model_tree.getListOfReactions();
        long number_of_reactions = listOfReactions.size();
        int enzyme_counter = 1;
        for (long reaction_index = 0;reaction_index<number_of_reactions;reaction_index++){

            // get the reaction object -
            Reaction reaction_object = listOfReactions.get(reaction_index);

            // What is the reaction name and Id?
            String reaction_id = reaction_object.getId();
            String reaction_name = reaction_object.getName();

            if (reaction_name.contains("Degradation") == false){

                // start writing -
                buffer.append("# ");
                buffer.append(reaction_name);
                buffer.append("\n");
                buffer.append("tmp_reaction = rate_constant_array[");
                buffer.append(reaction_index + 1);
                buffer.append("]*(E_");
                buffer.append(reaction_id);
                buffer.append(")");

                // Get the list of reactants -
                ListOfSpeciesReferences list_of_reactants = reaction_object.getListOfReactants();
                long number_of_reactants = list_of_reactants.size();
                for (long reactant_index = 0;reactant_index<number_of_reactants;reactant_index++){

                    // Get the species reference -
                    SimpleSpeciesReference species_reference = list_of_reactants.get(reactant_index);

                    // What is the id of this species?
                    String species_id = species_reference.getSpecies();
                    System.out.println("testing = "+species_id);

                    if (!species_id.equalsIgnoreCase("[]") && species_id.isEmpty() == false){

                        // What reactant index is this?
                        int local_reactant_index = tmp_species_symbol_vector.indexOf(species_id);

                        // Write the record -
                        buffer.append("*(");
                        buffer.append(species_id);
                        buffer.append("/(");
                        buffer.append("saturation_constant_array[");
                        buffer.append(reaction_index+1);
                        buffer.append(",");
                        buffer.append(local_reactant_index + 1);
                        buffer.append("] + ");
                        buffer.append(species_id);
                        buffer.append("))");
                    }
                }

                // add semicolon and new line -
                buffer.append(";\n");
                buffer.append("push!(rate_vector,tmp_reaction);\n");
                buffer.append("tmp_reaction = 0;\n");
                buffer.append("\n");
            }
            else {

                buffer.append("# ");
                buffer.append(reaction_name);
                buffer.append("\n");
                buffer.append("tmp_reaction = rate_constant_array[");
                buffer.append(reaction_index + 1);
                buffer.append("]*");
                buffer.append(reaction_id);

                // add semicolon and new line -
                buffer.append(";\n");
                buffer.append("push!(rate_vector,tmp_reaction);\n");
                buffer.append("tmp_reaction = 0;\n");
                buffer.append("\n");

                // update the enzyme counter -
                enzyme_counter++;
            }
        }

        // last line -
        buffer.append("# return the kinetics vector -\n");
        buffer.append("return rate_vector;\n");
        buffer.append("end\n");

        // return the buffer -
        return buffer.toString();
    }

    public String buildAllostericControlFunctionBuffer(Model model_tree, VLCGAllostericControlTreeWrapper control_tree,VLCGTransformationPropertyTree property_tree) throws Exception {

        // Method variables -
        StringBuffer buffer = new StringBuffer();

        // Get the control function name -
        String control_function_name = property_tree.lookupKwateeControlFunctionName();

        // Copyright notice -
        String copyright = copyrightFactory.getJuliaCopyrightHeader();
        buffer.append(copyright);

        // Fill in the buffer -
        buffer.append("function ");
        buffer.append(control_function_name);
        buffer.append("(t,x,rate_vector,data_dictionary)\n");
        buffer.append("# ---------------------------------------------------------------------- #\n");
        buffer.append("# ");
        buffer.append(control_function_name);
        buffer.append(".jl was generated using the Kwatee code generation system.\n");
        buffer.append("# Username: ");
        buffer.append(property_tree.lookupKwateeModelUsername());
        buffer.append("\n");
        buffer.append("# Type: ");
        buffer.append(property_tree.lookupKwateeModelType());
        buffer.append("\n");
        buffer.append("# Version: ");
        buffer.append(property_tree.lookupKwateeModelVersion());
        buffer.append("\n");
        buffer.append("# Generation timestamp: ");
        buffer.append(date_formatter.format(today));
        buffer.append("\n");
        buffer.append("# \n");
        buffer.append("# Arguments: \n");
        buffer.append("# t  - current time \n");
        buffer.append("# x  - state vector \n");
        buffer.append("# rate_vector - vector of reaction rates \n");
        buffer.append("# data_dictionary  - Data dictionary instance (holds model parameters) \n");
        buffer.append("# ---------------------------------------------------------------------- #\n");
        buffer.append("\n");
        buffer.append("# Set a default value for the allosteric control variables - \n");
        buffer.append("const number_of_reactions = length(rate_vector);\n");
        buffer.append("control_vector = ones(number_of_reactions);\n");
        buffer.append("const control_parameter_array = data_dictionary[\"CONTROL_PARAMETER_ARRAY\"];\n");
        buffer.append("\n");

        buffer.append("# Alias the species vector - \n");
        ListOfSpecies listOfSpecies = model_tree.getListOfSpecies();
        long number_of_species = listOfSpecies.size();
        for (long species_index = 0;species_index<number_of_species;species_index++){

            // Get the symbol -
            Species species_object = listOfSpecies.get(species_index);
            String species_symbol = species_object.getId();

            // write the symbol =
            //buffer.append("const ");
            buffer.append(species_symbol);
            buffer.append(" = x[");
            buffer.append(species_index + 1);
            buffer.append("];\n");
        }

        buffer.append("\n");

        // Use the control tree to formulate the control terms -
        // Process each reaction ... is this reaction a target?
        ListOfReactions reaction_list = model_tree.getListOfReactions();
        long number_of_reactions = reaction_list.size();
        int control_index = 1;
        for (long reaction_index = 0;reaction_index<number_of_reactions;reaction_index++){

            // Get the reaction object -
            Reaction reaction = reaction_list.get(reaction_index);
            String reaction_id = reaction.getId();

            // How many connections does this reaction have?
            Vector<VLCGAllostericControlModel> control_model_vector = control_tree.lookupAllostericConnectionsForReactionWithName(reaction_id);
            if (control_model_vector != null && !control_model_vector.isEmpty()){

                // ok, we have a non-empty vector of connections. We need to write terms for each
                Iterator<VLCGAllostericControlModel> control_model_iterator = control_model_vector.iterator();
                buffer.append("transfer_function_vector = Float64[];\n");
                while (control_model_iterator.hasNext()){

                    // Get the control model -
                    VLCGAllostericControlModel control_model_instance = control_model_iterator.next();

                    // I need to check the type -
                    String control_type = (String)control_model_instance.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_TYPE);
                    String control_actor = (String)control_model_instance.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_ACTOR);
                    String control_target = (String)control_model_instance.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_TARGET);

                    // Formulate the comment string -
                    String comment_string = "type: "+control_type+" actor: "+control_actor+" target: "+control_target;


                    buffer.append("# ");
                    buffer.append(comment_string);
                    buffer.append("\n");
                    if (control_type.equalsIgnoreCase("inhibition")){

                        // write -
                        buffer.append("push!(transfer_function_vector,1.0 - control_parameter_array[");
                        buffer.append(control_index);
                        buffer.append(",1]*");
                        buffer.append(control_model_instance.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_ACTOR));
                        buffer.append("^control_parameter_array[");
                        buffer.append(control_index);
                        buffer.append(",2]/(1+");
                        buffer.append("control_parameter_array[");
                        buffer.append(control_index);
                        buffer.append(",1]*");
                        buffer.append(control_model_instance.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_ACTOR));
                        buffer.append("^control_parameter_array[");
                        buffer.append(control_index);
                        buffer.append(",2]));\n");
                    }
                    else {

                        // write -
                        buffer.append("push!(transfer_function_vector,control_parameter_array[");
                        buffer.append(control_index);
                        buffer.append(",1]*");
                        buffer.append(control_model_instance.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_ACTOR));
                        buffer.append("^control_parameter_array[");
                        buffer.append(control_index);
                        buffer.append(",2]/(1+");
                        buffer.append("control_parameter_array[");
                        buffer.append(control_index);
                        buffer.append(",1]*");
                        buffer.append(control_model_instance.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_ACTOR));
                        buffer.append("^control_parameter_array[");
                        buffer.append(control_index);
                        buffer.append(",2]));\n");
                    }


                    // update index -
                    control_index++;
                }

                // ok, write the integration rule -
                buffer.append("control_vector[");
                buffer.append(reaction_index+1);
                buffer.append("] = mean(transfer_function_vector);\n");
                buffer.append("transfer_function_vector = 0;\n");
                buffer.append("\n");
            }
        }

        buffer.append("# Modify the rate_vector with the control variables - \n");
        buffer.append("rate_vector = rate_vector.*control_vector;\n");

        // last line -
        buffer.append("\n");
        buffer.append("# Return the modified rate vector - \n");
        buffer.append("return rate_vector;\n");
        buffer.append("end\n");

        // return the buffer -
        return buffer.toString();
    }

    public String buildBalanceFunctionBuffer(Model model_tree, VLCGTransformationPropertyTree property_tree) throws Exception {

        // Method variables -
        StringBuffer massbalances = new StringBuffer();

        // Get the balance function name -
        String balance_function_name = property_tree.lookupKwateeBalanceFunctionName();

        // Get/Set the kinetics function import -
        String kinetics_function_name = property_tree.lookupKwateeKineticsFunctionName();
        massbalances.append("include(\"");
        massbalances.append(kinetics_function_name);
        massbalances.append(".jl\");\n");

        // Get/Set the kinetics function import -
        String control_function_name = property_tree.lookupKwateeControlFunctionName();
        massbalances.append("include(\"");
        massbalances.append(control_function_name);
        massbalances.append(".jl\");\n");
        massbalances.append("\n");

        // Copyright notice -
        String copyright = copyrightFactory.getJuliaCopyrightHeader();
        massbalances.append(copyright);

        // Fill in the buffer -
        massbalances.append("function ");
        massbalances.append(balance_function_name);
        massbalances.append("(t,x,dxdt_vector,data_dictionary)\n");
        massbalances.append("# ---------------------------------------------------------------------- #\n");
        massbalances.append("# ");
        massbalances.append(balance_function_name);
        massbalances.append(".jl was generated using the Kwatee code generation system.\n");
        massbalances.append("# Username: ");
        massbalances.append(property_tree.lookupKwateeModelUsername());
        massbalances.append("\n");
        massbalances.append("# Type: ");
        massbalances.append(property_tree.lookupKwateeModelType());
        massbalances.append("\n");
        massbalances.append("# Version: ");
        massbalances.append(property_tree.lookupKwateeModelVersion());
        massbalances.append("\n");
        massbalances.append("# Generation timestamp: ");
        massbalances.append(date_formatter.format(today));
        massbalances.append("\n");
        massbalances.append("# \n");
        massbalances.append("# Arguments: \n");
        massbalances.append("# t  - current time \n");
        massbalances.append("# x  - state vector \n");
        massbalances.append("# dxdt_vector - right hand side vector \n");
        massbalances.append("# data_dictionary  - Data dictionary instance (holds model parameters) \n");
        massbalances.append("# ---------------------------------------------------------------------- #\n");
        massbalances.append("\n");
        massbalances.append("# Correct nagative x's = throws errors in control even if small - \n");
        massbalances.append("const idx = find(x->(x<0),x);\n");
        massbalances.append("x[idx] = 0.0;\n");
        massbalances.append("\n");
        massbalances.append("# Call the kinetics function - \n");
        massbalances.append("(rate_vector) = ");
        massbalances.append(kinetics_function_name);
        massbalances.append("(t,x,data_dictionary);\n");

        massbalances.append("\n");
        massbalances.append("# Call the control function - \n");
        massbalances.append("(rate_vector) = ");
        massbalances.append(control_function_name);
        massbalances.append("(t,x,rate_vector,data_dictionary);\n");
        massbalances.append("\n");

        // check - is this model large scale optimized?
        if (property_tree.isKwateeModelLargeScaleOptimized() == true){

            // build explicit list of balance equations -

        }
        else {

            // balance are encoded as matrix vector product -
            massbalances.append("# Encode the balance equations as a matrix vector product - \n");
            massbalances.append("const S = data_dictionary[\"STOICHIOMETRIC_MATRIX\"];\n");
            massbalances.append("const tmp_vector = S*rate_vector;\n");
            massbalances.append("const number_of_states = length(tmp_vector);\n");
            massbalances.append("for state_index in [1:number_of_states]\n");
            massbalances.append("\tdxdt_vector[state_index] = tmp_vector[state_index];\n");
            massbalances.append("end");
            massbalances.append("\n");
        }

        // last line -
        massbalances.append("\n");
        massbalances.append("end\n");

        // return the buffer -
        return massbalances.toString();
    }

    public String buildDriverFunctionBuffer(Model model_tree, VLCGTransformationPropertyTree property_tree) throws Exception {

        // String buffer -
        StringBuffer driver = new StringBuffer();

        // We need to get the imports -
        String balance_filename = property_tree.lookupKwateeBalanceFunctionName()+".jl";
        driver.append("include(\"");
        driver.append(balance_filename);
        driver.append("\")\n");
        driver.append("using Sundials;\n");
        driver.append("\n");

        // Copyright notice -
        String copyright = copyrightFactory.getJuliaCopyrightHeader();
        driver.append(copyright);

        // Get the function name -
        String function_name = property_tree.lookupKwateeDriverFunctionName();
        driver.append("function ");
        driver.append(function_name);
        driver.append("(TSTART,TSTOP,Ts,data_dictionary)\n");

        driver.append("# ----------------------------------------------------------------------------------- #\n");
        driver.append("# ");
        driver.append(function_name);
        driver.append(".jl was generated using the Kwatee code generation system.\n");
        driver.append("# ");
        driver.append(function_name);
        driver.append(": Solves model equations from TSTART to TSTOP given parameters in data_dictionary.\n");
        driver.append("# Username: ");
        driver.append(property_tree.lookupKwateeModelUsername());
        driver.append("\n");
        driver.append("# Type: ");
        driver.append(property_tree.lookupKwateeModelType());
        driver.append("\n");
        driver.append("# Version: ");
        driver.append(property_tree.lookupKwateeModelVersion());
        driver.append("\n");
        driver.append("# Generation timestamp: ");
        driver.append(date_formatter.format(today));
        driver.append("\n");
        driver.append("# \n");
        driver.append("# Input arguments: \n");
        driver.append("# TSTART  - Time start \n");
        driver.append("# TSTOP  - Time stop \n");
        driver.append("# Ts - Time step \n");
        driver.append("# data_dictionary  - Data dictionary instance (holds model parameters) \n");
        driver.append("# \n");
        driver.append("# Return arguments: \n");
        driver.append("# TSIM - Simulation time vector \n");
        driver.append("# X - Simulation state array (NTIME x NSPECIES) \n");
        driver.append("# ----------------------------------------------------------------------------------- #\n");
        driver.append("\n");

        driver.append("# Get required stuff from DataFile struct -\n");
        driver.append("TSIM = [TSTART:Ts:TSTOP];\n");
        driver.append("initial_condition_vector = data_dictionary[\"INITIAL_CONDITION_ARRAY\"];\n");
        driver.append("\n");

        driver.append("# Call the ODE solver - \n");
        driver.append("fbalances(t,y,ydot) = ");
        driver.append(property_tree.lookupKwateeBalanceFunctionName());
        driver.append("(t,y,ydot,data_dictionary);\n");
        driver.append("X = Sundials.cvode(fbalances,initial_condition_vector,TSIM);\n");
        driver.append("\n");
        driver.append("return (TSIM,X);\n");

        // last line -
        driver.append("end\n");

        // return the populated buffer -
        return driver.toString();
    }

    public String buildDataDictionaryBuffer(Model model_tree, VLCGAllostericControlTreeWrapper control_tree, VLCGTransformationPropertyTree property_tree) throws Exception {

        // String buffer -
        StringBuffer buffer = new StringBuffer();

        // Copyright notice -
        String copyright = copyrightFactory.getJuliaCopyrightHeader();
        buffer.append(copyright);

        // Get the function name -
        String function_name = property_tree.lookupKwateeDataDictionaryFunctionName();
        buffer.append("function ");
        buffer.append(function_name);
        buffer.append("(TSTART,TSTOP,Ts)\n");

        buffer.append("# ----------------------------------------------------------------------------------- #\n");
        buffer.append("# ");
        buffer.append(function_name);
        buffer.append(".jl was generated using the Kwatee code generation system.\n");
        buffer.append("# ");
        buffer.append(function_name);
        buffer.append(": Stores model parameters as key - value pairs in a Julia Dict() \n");
        buffer.append("# Username: ");
        buffer.append(property_tree.lookupKwateeModelUsername());
        buffer.append("\n");
        buffer.append("# Type: ");
        buffer.append(property_tree.lookupKwateeModelType());
        buffer.append("\n");
        buffer.append("# Version: ");
        buffer.append(property_tree.lookupKwateeModelVersion());
        buffer.append("\n");
        buffer.append("# Generation timestamp: ");
        buffer.append(date_formatter.format(today));
        buffer.append("\n");
        buffer.append("# \n");
        buffer.append("# Input arguments: \n");
        buffer.append("# TSTART  - Time start \n");
        buffer.append("# TSTOP  - Time stop \n");
        buffer.append("# Ts - Time step \n");
        buffer.append("# \n");
        buffer.append("# Return arguments: \n");
        buffer.append("# data_dictionary  - Data dictionary instance (holds model parameters) \n");
        buffer.append("# ----------------------------------------------------------------------------------- #\n");
        buffer.append("\n");
        buffer.append("# Load the stoichiometric matrix - \n");

        // Get the path to the stoichiometric matrix -
        String fully_qualified_stoichiometric_matrix_path = property_tree.lookupKwateeStoichiometricMatrixFilePath();
        buffer.append("S = float(open(readdlm,");
        buffer.append("\"");
        buffer.append(fully_qualified_stoichiometric_matrix_path);
        buffer.append("\"));\n");
        buffer.append("(NSPECIES,NREACTIONS) = size(S);\n");

        buffer.append("\n");
        buffer.append("# Formulate the initial condition array - \n");
        buffer.append("initial_condition_array = Float64[];\n");
        ListOfSpecies listOfSpecies = model_tree.getListOfSpecies();
        Vector<String> tmp_species_symbol_vector = new Vector<String>();
        long number_of_species = listOfSpecies.size();
        for (long species_index = 0;species_index<number_of_species;species_index++){

            // Get the species at index species_index -
            Species species_object = listOfSpecies.get(species_index);
            String species_symbol = species_object.getId();

            // What is the ic for this species?
            double ic_value = species_object.getInitialConcentration();

            // write ic record -
            buffer.append("push!(initial_condition_array,");
            buffer.append(ic_value);
            buffer.append(");\t");
            buffer.append("#\t");
            buffer.append(species_index+1);
            buffer.append("\t id:");
            buffer.append(species_object.getId());
            buffer.append("\t symbol:");
            buffer.append(species_object.getName());
            buffer.append("\n");

            // store the species id -
            tmp_species_symbol_vector.addElement(species_symbol);
        }

        buffer.append("\n");
        buffer.append("# Formulate the rate constant array - \n");
        buffer.append("rate_constant_array = Float64[];\n");
        ListOfReactions listOfReactions = model_tree.getListOfReactions();
        long number_of_reactions = listOfReactions.size();
        for (long reaction_index = 0;reaction_index<number_of_reactions;reaction_index++){

            // Get the reaction object at reaction_index -
            Reaction reaction_object = listOfReactions.get(reaction_index);

            // What is the name of this reaction?
            String reaction_name = reaction_object.getName();

            // write parameter record -
            buffer.append("push!(rate_constant_array,0.0)\t");
            buffer.append("#\t");
            buffer.append(reaction_index + 1);
            buffer.append("\t ");
            buffer.append(reaction_name);
            buffer.append("\n");
        }

        // Saturation constants -
        buffer.append("\n");
        buffer.append("# Formulate the saturation constant array - \n");
        buffer.append("saturation_constant_array = zeros(NREACTIONS,NSPECIES);\n");
        for (long reaction_index = 0;reaction_index<number_of_reactions;reaction_index++){

            // Get the reaction object at reaction_index -
            Reaction reaction_object = listOfReactions.get(reaction_index);

            // Get the reaction name/id -
            String reaction_id = reaction_object.getName();

            ListOfSpeciesReferences list_of_reactants = reaction_object.getListOfReactants();
            long number_of_reactants = list_of_reactants.size();
            for (long reactant_index = 0;reactant_index<number_of_reactants;reactant_index++) {

                // Get the species reference -
                SimpleSpeciesReference species_reference = list_of_reactants.get(reactant_index);

                // What is the id of this species?
                String species_id = species_reference.getSpecies();
                if (!species_id.equalsIgnoreCase("[]") &&
                        species_id.isEmpty() == false &&
                        !reaction_id.contains("Degradation")) {

                    // What reactant index is this?
                    int local_reactant_index = tmp_species_symbol_vector.indexOf(species_id);

                    // ok, we should have all we need to replace the non-zero values for the K array -
                    buffer.append("saturation_constant_array[");
                    buffer.append(reaction_index+1);
                    buffer.append(",");
                    buffer.append(local_reactant_index+1);
                    buffer.append("] = 1.0;");
                    buffer.append("\t#\t Name: ");
                    buffer.append(reaction_id);
                    buffer.append(" Species: ");
                    buffer.append(species_id);
                    buffer.append("\n");
                }
            }
        }

        // Control model constants -
        buffer.append("\n");
        buffer.append("# Formulate control parameter array - \n");
        int number_of_control_terms = control_tree.calculateNumberOfControlTerms();
        buffer.append("control_parameter_array = zeros(");
        buffer.append(number_of_control_terms);
        buffer.append(",2);\n");

        int control_term_index = 1;
        for (long reaction_index = 0;reaction_index<number_of_reactions;reaction_index++) {

            // Get the reaction object -
            Reaction reaction = listOfReactions.get(reaction_index);
            String reaction_id = reaction.getId();

            // How many connections does this reaction have?
            Vector<VLCGAllostericControlModel> control_model_vector = control_tree.lookupAllostericConnectionsForReactionWithName(reaction_id);
            if (control_model_vector != null && !control_model_vector.isEmpty()) {

                Iterator<VLCGAllostericControlModel> control_model_iterator = control_model_vector.iterator();
                while (control_model_iterator.hasNext()) {

                    // Get the control model -
                    VLCGAllostericControlModel controlModel = control_model_iterator.next();

                    // Get items from the model required for the comment -
                    String control_type = (String)controlModel.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_TYPE);
                    String control_actor = (String)controlModel.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_ACTOR);
                    String control_target = (String)controlModel.getAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_TARGET);

                    // Formulate the comment string -
                    String comment_string = "type: "+control_type+" actor: "+control_actor+" target: "+control_target;

                    // write the entries in the control array -
                    buffer.append("control_parameter_array["+(control_term_index)+",1] = 0.1;\t # gain => "+comment_string+"\n");
                    buffer.append("control_parameter_array["+(control_term_index)+",2] = 1.0;\t # order => "+comment_string+"\n");

                    // update the counter -
                    control_term_index++;
                }
            }
        }


        buffer.append("\n");
        buffer.append("# ---------------------------- DO NOT EDIT BELOW THIS LINE -------------------------- #\n");
        buffer.append("data_dictionary = Dict();\n");
        buffer.append("data_dictionary[\"STOICHIOMETRIC_MATRIX\"] = S;\n");
        buffer.append("data_dictionary[\"RATE_CONSTANT_ARRAY\"] = rate_constant_array;\n");
        buffer.append("data_dictionary[\"SATURATION_CONSTANT_ARRAY\"] = saturation_constant_array;\n");
        buffer.append("data_dictionary[\"INITIAL_CONDITION_ARRAY\"] = initial_condition_array;\n");
        buffer.append("data_dictionary[\"CONTROL_PARAMETER_ARRAY\"] = control_parameter_array;\n");
        buffer.append("# ----------------------------------------------------------------------------------- #\n");

        // last line -
        buffer.append("return data_dictionary;\n");
        buffer.append("end\n");

        // return the populated buffer -
        return buffer.toString();
    }

    public void buildStoichiometricMatrix(double[][] dblSTMatrix,Model model_wrapper) throws Exception {

        // Get the dimension of the system -
        int NUMBER_OF_SPECIES = 0;
        int NUMBER_OF_RATES = 0;

        // Get the system dimension -
        NUMBER_OF_SPECIES = (int)model_wrapper.getNumSpecies();
        NUMBER_OF_RATES = (int)model_wrapper.getNumReactions();

        // Go through and put everything as zeros by default -
        for (int scounter=0;scounter<NUMBER_OF_SPECIES;scounter++)
        {
            for (int rcounter=0;rcounter<NUMBER_OF_RATES;rcounter++)
            {
                dblSTMatrix[scounter][rcounter]=0.0;
            }
        }

        // When I get here, I have a st. matrix w/all zeros -
        // put in the correct values -
        ListOf listRates = model_wrapper.getListOfReactions();
        ListOf listSpecies = model_wrapper.getListOfSpecies();

        for (int scounter=0;scounter<NUMBER_OF_SPECIES;scounter++)
        {
            // Get the species reference -
            Species species = (Species)listSpecies.get(scounter);
            String strSpecies = species.getName();

            // Ok, I need to go through the rates and determine if this species is involved -
            for (int rcounter=0;rcounter<NUMBER_OF_RATES;rcounter++)
            {
                // Get the Reaction object -
                Reaction rxn_obj = (Reaction)listRates.get(rcounter);

                // Get the 'radius' of this rate -
                int NUMBER_OF_REACTANTS = (int)rxn_obj.getNumReactants();
                int NUMBER_OF_PRODUCTS = (int)rxn_obj.getNumProducts();

                // go through the reactants of this reaction -
                for (int reactant_index=0;reactant_index<NUMBER_OF_REACTANTS;reactant_index++)
                {
                    // Get the species reference -
                    SpeciesReference species_ref = rxn_obj.getReactant(reactant_index);
                    String strReactant = species_ref.getSpecies();

                    if (strReactant.equalsIgnoreCase(strSpecies))
                    {

                        dblSTMatrix[scounter][rcounter]=-1.0*species_ref.getStoichiometry();
                    }
                }

                // go through the products of this reaction -
                for (int product_index=0;product_index<NUMBER_OF_PRODUCTS;product_index++)
                {
                    // Get the species reference -
                    SpeciesReference species_ref = rxn_obj.getProduct(product_index);
                    String strProduct = species_ref.getSpecies();

                    //System.out.println("Comparing NP="+NUMBER_OF_PRODUCTS+" to "+strProduct+"="+strSpecies+"?");

                    if (strProduct.equalsIgnoreCase(strSpecies))
                    {
                        dblSTMatrix[scounter][rcounter]=species_ref.getStoichiometry();
                    }
                }
            }
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
