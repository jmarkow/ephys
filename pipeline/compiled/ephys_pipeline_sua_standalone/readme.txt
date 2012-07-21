MATLAB Compiler

1. Prerequisites for Deployment 

. Verify the MATLAB Compiler Runtime (MCR) is installed and ensure you    
  have installed version 7.16.   

. If the MCR is not installed, launch MCRInstaller, located in:

  <matlabroot>*/toolbox/compiler/deploy/maci64/MCRInstaller.zip

For more information about the MCR and the MCR Installer, see 
“Working With the MCR” in the MATLAB Compiler User’s Guide.    


NOTE: You will need administrator rights to run MCRInstaller. 


2. Files to Deploy and Package

Files to package for Standalone 
================================
-run_ephys_pipeline_sua_standalone.sh (shell script run to temporarily set environment 
 variables and execute the application)
   -to run the shell script, type
   
       ./run_ephys_pipeline_sua_standalone.sh <mcr_directory> <argument_list>
       
    at Linux or Mac command prompt. <mcr_directory> is the directory 
    where version 7.16 of MCR is installed or the directory where 
    MATLAB is installed on the machine. <argument_list> is all the 
    arguments you want to pass to your application. For example, 

    If you have version 7.16 of MCR installed in 
    /mathworks/home/application/R2010a/v716, run the shell script as:
    
       ./run_ephys_pipeline_sua_standalone.sh /mathworks/home/application/R2010a/v716
       
    If you have MATLAB installed in /mathworks/devel/application/matlab, 
    run the shell script as:
    
       ./run_ephys_pipeline_sua_standalone.sh /mathworks/devel/application/matlab
-MCRInstaller.zip 
   -include when building component by clicking "Add MCR" link 
    in deploytool
-The Macintosh bundle directory structure ephys_pipeline_sua_standalone.app 
   -this can be gathered up using the zip command 
    zip -r ephys_pipeline_sua_standalone.zip ephys_pipeline_sua_standalone.app
    or the tar command 
    tar -cvf ephys_pipeline_sua_standalone.tar ephys_pipeline_sua_standalone.app
-This readme file 

3. Definitions

For information on deployment terminology, go to 
http://www.mathworks.com/help. Select your product and see 
the Glossary in the User’s Guide.


* NOTE: <matlabroot> is the directory where MATLAB is installed on the target machine.


4. Appendix 

A. Mac systems:
   On the target machine, add the MCR directory to the environment variable 
   DYLD_LIBRARY_PATH by issuing the following commands:

        NOTE: <mcr_root> is the directory where MCR is installed
              on the target machine.         

            setenv DYLD_LIBRARY_PATH
                $DYLD_LIBRARY_PATH:
                <mcr_root>/v716/runtime/maci64:
                <mcr_root>/v716/sys/os/maci64:
                <mcr_root>/v716/bin/maci64:
                /System/Library/Frameworks/JavaVM.framework/JavaVM:
                /System/Library/Frameworks/JavaVM.framework/Libraries
            setenv XAPPLRESDIR <mcr_root>/v716/X11/app-defaults



     
        NOTE: To make these changes persistent after logout on Linux 
              or Mac machines, modify the .cshrc file to include this  
              setenv command.
        NOTE: The environment variable syntax utilizes forward 
              slashes (/), delimited by colons (:).  
        NOTE: When deploying standalone applications, it is possible 
              to run the shell script file run_ephys_pipeline_sua_standalone.sh 
              instead of setting environment variables. See 
              section 2 "Files to Deploy and Package".    



5. Launching of application using Macintosh finder.

If the application is purely graphical, that is, it doesn't read from standard in or 
write to standard out or standard error, it may be launched in the finder just like any 
other Macintosh application.



