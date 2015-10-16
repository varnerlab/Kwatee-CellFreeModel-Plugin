package org.varnerlab.kwatee.cellfreemodel;

import org.w3c.dom.Document;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathFactory;
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
 * Created by jeffreyvarner on 10/16/15.
 */
public class VLCGAllostericControlTreeWrapper {

    // instance variables -
    private Document _allosteric_control_tree = null;
    private XPathFactory _xpath_factory = XPathFactory.newInstance();
    private XPath _xpath = _xpath_factory.newXPath();

    public VLCGAllostericControlTreeWrapper() {
    }

    public void set_allosteric_control_tree(Document allosteric_control_tree) {
        this._allosteric_control_tree = allosteric_control_tree;
    }



    public int calculateNumberOfControlTerms() throws Exception {

        // method variables -
        int number_of_control_terms = 0;

        // get the vector of connections -
        String xpath_string = ".//allosteric/connection";

        // Lookup -
        NodeList connectionNodeList = (NodeList) _xpath.evaluate(xpath_string, _allosteric_control_tree, XPathConstants.NODESET);
        number_of_control_terms = connectionNodeList.getLength();

        // return -
        return number_of_control_terms;
    }

    public Vector<VLCGAllostericControlModel> lookupAllostericConnectionsForReactionWithName(String reaction_name) throws Exception {

        if (reaction_name == null || _allosteric_control_tree == null){
            throw new Exception("ERROR: No reaction name was requested, or we are missing the allosteric control tree?");
        }


        // method variables -
        Vector<VLCGAllostericControlModel> model_vector = null;

        // Formulate the xpath_string -
        String xpath_string = ".//allosteric/connection[@target='"+reaction_name+"']";

        // Lookup -
        NodeList connectionNodeList = (NodeList) _xpath.evaluate(xpath_string, _allosteric_control_tree, XPathConstants.NODESET);
        int number_of_connections = connectionNodeList.getLength();

        // do we have connections?
        if (number_of_connections>0){
            model_vector = new Vector<VLCGAllostericControlModel>();
        }

        for (int connection_index = 0;connection_index<number_of_connections;connection_index++){

            // Get a node ..
            Node connection_node = connectionNodeList.item(connection_index);

            // Get the attributes for this nodel -
            NamedNodeMap attribute_map = connection_node.getAttributes();
            Node type_node = attribute_map.getNamedItem("type");
            Node actor_node = attribute_map.getNamedItem("actor");
            Node name_node = attribute_map.getNamedItem("name");

            // Create a new control model -
            VLCGAllostericControlModel control_model = new VLCGAllostericControlModel();
            control_model.setAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_NAME,name_node.getNodeValue());
            control_model.setAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_TARGET,reaction_name);
            control_model.setAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_TYPE,type_node.getNodeValue());
            control_model.setAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_ACTOR,actor_node.getNodeValue());
            model_vector.addElement(control_model);
        }

        // return -
        return model_vector;
    }

}
