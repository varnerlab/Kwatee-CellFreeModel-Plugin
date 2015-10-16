package org.varnerlab.kwatee.cellfreemodel;

// imports -
import java.util.Hashtable;
import java.util.Vector;
import java.util.StringTokenizer;
import java.util.Iterator;

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

public class VLCGMetabolicReactionModel {

    // Instance variables -
    private Hashtable _reaction_component_table = new Hashtable();
    private Vector _vecReactants;
    private Vector _vecProducts;
    private String _strReverseFlag = "";
    private String _strForwardFlag = "";
    private String _strFormatedRxnString = "";

    // keys -
    public static final String REACTION_NAME = "REACTION_NAME";
    public static final String REACTION_REACTANTS = "REACTANTS";
    public static final String REACTION_PRODUCTS = "PRODUCTS";
    public static final String REACTION_REVERSE = "REVERSE";
    public static final String REACTION_FORWARD = "FORWARD";
    public static final String REACTION_TYPE = "RTYPE";

    public VLCGMetabolicReactionModel() {
        _init();
    }

    private void _init() {

        // Create new instances of the reactants and products
        _vecReactants = new Vector();
        _vecProducts = new Vector();
    }

    public void setMetabolicReactionComponent(String key,Object value) {

        if (key == null && value == null){
            return;
        }

        System.out.println("Storing key = "+key+" value = "+value);

        // store the reaction component -
        _reaction_component_table.put(key,value);
    }

    public String toString(){

        String debug_string = "";
        return _strFormatedRxnString+" R:"+_strReverseFlag+" F:"+_strForwardFlag;

    }

    public Object getMetabolicReactionComponent(String key) throws Exception {

        if (key == null || _reaction_component_table.containsKey(key) == false){
            throw new Exception("ERROR: Missing metabolic reaction component. Can't find key = "+key);
        }

        return _reaction_component_table.get(key);
    }

    public void doExecute() throws Exception {

        // Get the reactants and products -
        String reaction_name = (String)getMetabolicReactionComponent(REACTION_NAME);
        String reactant_string = (String)getMetabolicReactionComponent(REACTION_REACTANTS);
        String product_string = (String)getMetabolicReactionComponent(REACTION_PRODUCTS);
        _doExecute(reactant_string,product_string);

        // formulate the string -
        _strFormatedRxnString = reaction_name+": "+reactant_string+" = "+product_string;
    }

    private void _doExecute(String rxnStringReactants,String rxnStringProducts) throws Exception {

        // Ok, so I have the reaction strings, cut them apart
        parseString(rxnStringReactants,_vecReactants,false);
        parseString(rxnStringProducts,_vecProducts,true);
    }

    public double getCoefficientValue(String symbol) {

        // Method attributes
        double coeff = 0.0;

        // Ok, so here is the hard part -

        // Go through the reactants -
        Iterator iterReactants = this.getReactants();
        while (iterReactants.hasNext()) {
            // Get the state symbol
            VLCGMetaboliteSymbol tmp = (VLCGMetaboliteSymbol)iterReactants.next();
            String tmpString = tmp.getSymbol();

            if (tmpString.equalsIgnoreCase(symbol))
            {
                // If I get here then I have a match -
                coeff = tmp.getCoefficient();
            }
        }

        // Go through the products -
        Iterator iterProducts = this.getProducts();
        while (iterProducts.hasNext())
        {
            // Get the state symbol
            VLCGMetaboliteSymbol tmp = (VLCGMetaboliteSymbol)iterProducts.next();
            String tmpString = tmp.getSymbol();

            if (tmpString.equalsIgnoreCase(symbol))
            {
                // If I get here then I have a match -
                coeff = tmp.getCoefficient();
            }
        }

        // return coeff to the caller
        return(coeff);
    }

    public String getBoundsString() {
        String strTmp = "";

        strTmp = this._strReverseFlag+"\t"+this._strForwardFlag;

        return(strTmp);
    }

    public String getReactionString() {
        return(this._strFormatedRxnString);
    }

    // Get methods for the reactants and products vectors
    public Iterator getReactants() {
        return(_vecReactants.iterator());
    }

    public Iterator getProducts() {
        return(_vecProducts.iterator());
    }

    private void parseString(String frag,Vector vector,boolean isProduct) throws Exception {
        // Ok, this method contains the logic to cut up the reaction strings -

        // Cut around the +'s'
        StringTokenizer tokenizer=new StringTokenizer(frag,"+",false);
        while (tokenizer.hasMoreElements())
        {
            // Get a data from the tokenizer -
            Object dataChunk=tokenizer.nextToken();

            // Create new symbol wrapper
            VLCGMetaboliteSymbol symbol = new VLCGMetaboliteSymbol();

            // Check to see if this dataChunk string contains a *
            if (((String)dataChunk).contains("*"))
            {
                // If I get here, then the string contains a stoichiometric coefficient

                // Cut around the *'s
                StringTokenizer tokenizerCoeff=new StringTokenizer((String)dataChunk,"*",false);
                int intCoeffCounter = 1;
                while (tokenizerCoeff.hasMoreElements())
                {
                    Object dataCoeff=tokenizerCoeff.nextToken();

                    if (intCoeffCounter==1)
                    {
                        if (isProduct){
                            symbol.setCoefficient(Double.parseDouble(((String)dataCoeff)));
                        }
                        else
                        {
                            double dblTmp = Double.parseDouble(((String)dataCoeff));
                            symbol.setCoefficient(-1.0*dblTmp);
                        }

                        // Update the counter
                        intCoeffCounter++;
                    }
                    else if (intCoeffCounter==2)
                    {
                        symbol.setSymbol((String)dataCoeff);
                    }
                }
            }
            else
            {
                // If I get here, then no coefficient
                if (isProduct)
                {
                    // If this metabolite is in a product string, then coeff is positive
                    symbol.setSymbol((String)dataChunk);
                    symbol.setCoefficient(1.0);
                }
                else
                {
                    // If this metabolite is in a reactant string, then coeff is negative
                    symbol.setSymbol((String)dataChunk);
                    symbol.setCoefficient(-1.0);
                }

            }

            // Add to symbol wrapper to the vector -
            vector.addElement(symbol);
        }
    }
}
