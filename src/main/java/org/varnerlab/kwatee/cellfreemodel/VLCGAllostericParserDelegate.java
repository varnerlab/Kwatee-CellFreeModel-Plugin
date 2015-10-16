package org.varnerlab.kwatee.cellfreemodel;

import java.io.StringReader;
import java.util.StringTokenizer;

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
 * Created by jeffreyvarner on 10/12/15.
 */

public class VLCGAllostericParserDelegate implements VLCGParserHandlerDelegate {

    public VLCGAllostericParserDelegate() {
    }


    @Override
    public Object parseLine(String line) throws Exception {

        // Build the control model -
        VLCGAllostericControlModel control_model = new VLCGAllostericControlModel();

        StringTokenizer tokenizer=new StringTokenizer(cleanRecord(line),",",false);
        int intCounter = 1;
        while (tokenizer.hasMoreElements()) {

            // Get a data from the tokenizer -
            Object dataChunk = tokenizer.nextToken();

            if (intCounter == 1) {
                control_model.setAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_NAME, dataChunk);
            }
            else if (intCounter == 2){
                control_model.setAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_ACTOR, dataChunk);
            }
            else if (intCounter == 3){
                control_model.setAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_TARGET, dataChunk);
            }
            else if (intCounter == 4){

                // cut off the trailing ;
                String type_string = (String)dataChunk;
                String cut_string = type_string.substring(0,type_string.length()-1);

                control_model.setAllostericControlComponent(VLCGAllostericControlModel.ALLOSTERIC_CONTROL_TYPE, cut_string);
            }

            // update the counter -
            intCounter++;
        }

        return control_model;
    }

    private String cleanRecord(String record) {
        String rString = record.replaceAll(",,",",NULL,");
        return(rString);
    }
}
