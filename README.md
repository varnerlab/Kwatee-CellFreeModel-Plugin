Introduction
----

The Kwatee cell free model plugin generates a dynamic model of cell free metabolism using the hybrid modeling framework of Wayman et al:

[Wayman J, Sagar A and J. Varner (2015) Dynamic Modeling of Cell Free Biochemical Networks using Effective Kinetic Models. Processes DOI:10.3390/pr3010138](http://www.mdpi.com/2227-9717/3/1/138)

The cell free model equations are generated in the [Julia](http://julialang.org) programming language. The model equations are solved using the [SUNDIALS](https://github.com/JuliaLang/Sundials.jl/blob/master/README.md) package for [Julia](http://julialang.org). 

How do I generate code?
---

The plugin can be used by putting the jar file into the plugins subdirectory of your [Kwatee server installation](https://github.com/varnerlab/KwateeServer). You can either use the current compiled jar file found at: 

~~~
./build/libs
~~~

or you can compile the source code yourself using the [Gradle](http://gradle.org) build system by executing:

~~~
gradle jar
~~~

in the root directory of this repository. __Note:__ If you choose to compile from source, you need both the Kwatee server and libSBML jars on your Java classpath. The `build.gradle` file in the repository automatically adds all jars found in the `libs` subdirectory to your classpath.


Model and job files
----

The biology of your cell free model is specified using a simple comma delimated file structure called a Varner flat file (VFF). An example VFF for a cell free model is included in the `examples` subdirectory of the [Kwatee server](https://github.com/varnerlab/KwateeServer/tree/master/examples/cell-free-example) and shown below:

~~~
// ===================================================================== //
// Test model file for the Kwatee cell free plugin
// Author: J. Varner
// School of Chemical Engineering, Purdue University
// West Lafayette, IN 47097
//
// Version: 1.01
// ==================================================================== //

// ==================================================================== //
// Metabolic reactions -
[kwatee.metabolic.reaction.handler]
//
// Reaction records are structured as:
// NAME (unique),reactants,products,reverse {0|-inf},forward {0|inf};
// ==================================================================== //

R_1,A+2*B,C,0,inf;
R_2,C,D+E,0,inf;
R_3,D,[],0,inf;
R_4,E,F,-inf,inf;
R_5,F,[],0,inf;
R_6,[],A,0,inf;
R_7,[],B,0,inf;

// ==================================================================== //
// Allosteric regulation -
[kwatee.allosteric.control.handler]
//
// Allosteric records have the form:
// NAME (unique),actor (regulator),target,type {inhibition|activation}
// ==================================================================== //
INHIBITION_RXN_NAME,C,R_6,inhibition;
ACTIVATION_RXN_NAME_2,A,R_1,activation;
ACTIVATION_RXN_NAME,A,R_4,activation;
~~~

The model specification file (by default given the filename `Model.net`) defines the biology of the model that will get generated. A cell free VFF contains two sections, the metabolic reaction section (top) and the allosteric regulation section (bottom). In addition to the model specification, to succesfully generate code you will need a job specification file by default given the filename `Configuration.xml`. A typical job configuration file is:

~~~
<?xml version="1.0" encoding="UTF-8"?>
<Model username="jeffreyvarner" model_version="1.0" model_type="CFPS-JULIA" large_scale_optimized="false" model_name="TEST_MODEL">
  <Configuration>
    <ListOfPackages>
        <package required="YES" symbol="INPUT_HANDLER_PACKAGE" package_name="org.varnerlab.kwatee.cellfreemodel"></package>
        <package required="YES" symbol="OUTPUT_HANDLER_PACKAGE" package_name="org.varnerlab.kwatee.cellfreemodel"></package>
    </ListOfPackages>
    <ListOfPaths>
        <path required="YES" symbol="KWATEE_INPUT_PATH" path_location="/Users/jeffreyvarner/Desktop/development/kwatee/examples/cell-free-example/"></path>
        <path required="YES" symbol="KWATEE_SOURCE_OUTPUT_PATH" path_location="/Users/jeffreyvarner/Desktop/development/kwatee/examples/cell-free-example/src/"></path>
        <path required="YES" symbol="KWATEE_NETWORK_OUTPUT_PATH" path_location="/Users/jeffreyvarner/Desktop/development/kwatee/examples/cell-free-example/network/"></path>
        <path required="YES" symbol="KWATEE_DEBUG_OUTPUT_PATH" path_location="/Users/jeffreyvarner/Desktop/development/kwatee/examples/cell-free-example/debug/"></path>

        <path required="YES" symbol="KWATEE_SERVER_ROOT_DIRECTORY" path_location="/Users/jeffreyvarner/Desktop/KWATEEServer-v1.0/"></path>
        <path required="YES" symbol="KWATEE_SERVER_JAR_DIRECTORY" path_location="/Users/jeffreyvarner/Desktop/KWATEEServer-v1.0/dist/"></path>
        <path required="YES" symbol="KWATEE_PLUGINS_JAR_DIRECTORY" path_location="/Users/jeffreyvarner/Desktop/KWATEEServer-v1.0/plugins/"></path>
    </ListOfPaths>
  </Configuration>

  <Handler>
      <InputHandler required="YES" input_classname="VLCGParseVarnerFlatFile" package="INPUT_HANDLER_PACKAGE"></InputHandler>
      <OutputHandler required="YES" output_classname="VLCGWriteJuliaCellFreeModel" package="OUTPUT_HANDLER_PACKAGE"></OutputHandler>
  </Handler>
  <InputOptions>
      <NetworkFile required="YES" path_symbol="KWATEE_INPUT_PATH" filename="Model.net"></NetworkFile>
      <OrderFile required="NO" path_symbol="KWATEE_INPUT_PATH" filename="Order.dat"></OrderFile>
      <ModelParameterFile required="NO" path_symbol="KWATEE_INPUT_PATH" filename="Parameters.dat"></ModelParameterFile>
      <InitialConditionFile required="NO" path_symbol="KWATEE_INPUT_PATH" filename="InitialConditins.dat"></InitialConditionFile>
  </InputOptions>
  <OutputOptions>
      <DataFunction required="YES" path_symbol="KWATEE_SOURCE_OUTPUT_PATH" filename="DataFile.jl"></DataFunction>
      <BalanceFunction required="YES" path_symbol="KWATEE_SOURCE_OUTPUT_PATH" filename="MassBalances.jl"></BalanceFunction>
      <KineticsFunction required="YES" path_symbol="KWATEE_SOURCE_OUTPUT_PATH" filename="Kinetics.jl"></KineticsFunction>
      <InputFunction required="YES" path_symbol="KWATEE_SOURCE_OUTPUT_PATH" filename="Inputs.jl"></InputFunction>
      <DriverFunction required="YES" path_symbol="KWATEE_SOURCE_OUTPUT_PATH" filename="SolveBalances.jl"></DriverFunction>
      <ControlFunction required="YES" path_symbol="KWATEE_SOURCE_OUTPUT_PATH" filename="Control.jl"></ControlFunction>

      <StoichiometricMatrix required="YES" path_symbol="KWATEE_NETWORK_OUTPUT_PATH" filename="Network.dat"></StoichiometricMatrix>
      <DebugOutputFile required="YES" path_symbol="KWATEE_DEDUG_OUTPUT_PATH" filename="Debug.txt"></DebugOutputFile>
  </OutputOptions>
</Model>
~~~

The majority of the fields in the job configuration file can stay at the default values. However, you will need to specify your path structure (where Kwatee can find your files, where you want your generated code to reside, and where to find the server). Thus, you should edit the paths in the `<listOfPaths>...</listOfPaths>` section of the configuration file with your values. There typically only a few paths that must be specified:

* `KWATEE_INPUT_PATH`: Directory where Kwatee will find your `Model.net` file.
* `KWATEE_SOURCE_OUTPUT_PATH`: Directory where your generated model source code will be written (default is `src`).
* `KWATEE_NETWORK_OUTPUT_PATH`: Directory where your stoichiometric matrix be written (default is `network`).
* `KWATEE_DEBUG_OUTPUT_PATH`: Directory where any debug information is written (default is `debug`).
* `KWATEE_SERVER_ROOT_DIRECTORY`: Directory where your Kwatee server is installed.
* `KWATEE_SERVER_JAR_DIRECTORY`: Subdirectory where the Kwatee server jar file can be found (default is `dist`).
* `KWATEE_PLUGINS_JAR_DIRECTORY`: Subdirectory where Kwatee can find your plugin jars (default is `plugins`).

