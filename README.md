### Container for under-clustering detection and correction of gene-famiiies

This container contains the application for detecting and correcting under-clustered gene families using a sequence-pair-based classification method. Briefly, the method first trains and tests HMM models based on sequence pairs from a given family fasta and the candidate missing sequences(a. k. a. closest non-family sequences) for the family, to determine the separation between the family sequences and the non-family sequences. Classification statistics between the family and non-family sequences are then used to predict missing sequences for the family from the list of candidate missing sequences.

#### Steps for executing the analysis
This application can be executed in 2 modes: global and subset. The global mode is used to predict missing sequences for a given family by searching in all the rest of the families. On the other hand, the subset mode is used to predict missing sequences from a non-overlapping subset of families (typically usually small families that could be merged into larger families). Further, both types of modes can be executed with 2 prediction modes: F1-score and F2-score. The F1-score mode is the conservative mode for prediction of missing sequences that places equal weights on precision and recall while predicting missing sequences for families. The F2-score mode favors recall over precision for predicting missing sequences and thus can predict more missing sequences than the F1-score mode for the same family. If the user is confident that the given families are highly accurate with low probability of them containing "wrong" or unrelated sequences, F2-score mode can be used for correcting the families. Else, if the user is not confident about the accuracy of the families, the F1-score mode is the safest to avoid attracting more unrelated sequences to the families.

 1. #### Downloading the container
  ```
  docker pull akshayayadav/overcl-detection-correction
  ```

 2. #### Preparing the data
  Prepare a data directory *<my_directory>* with a user defined name containing a directory named *family_fasta* and a fasta file from which to search and predict missing sequences for the families. Fasta files for all the families to be analyzed should be placed in the *family_fasta* directory. For global mode, the fasta file should be named *proteomes.fa* and **MUST** contain sequences from all the families including the families present in the *family_fasta* directory. For subset mode, the fasta file should be named *unclustered.fa* and **MUST NOT** contain any sequences from the families present in *family_fasta* directory.

 3. #### Running the analysis
  * Global mode with F1-score function and *\<n\>* cores
  ```
  docker run -v <absolute_path_to_data_directory>:/data akshayayadav/overcl-detection-correction run_analysis_global-F1.sh -c <n>
  ```
  
  * Global mode with F2-score function and *\<n\>* cores
  ```
  docker run -v <absolute_path_to_data_directory>:/data akshayayadav/overcl-detection-correction run_analysis_global-F2.sh -c <n>
  ```
  
  * Subset mode with F1-score function and *\<n\>* cores
  ```
  docker run -v <absolute_path_to_data_directory>:/data akshayayadav/overcl-detection-correction run_analysis_subset-F1.sh -c <n>
  ```
  
  * Subset mode with F2-score function and *\<n\>* cores
  ```
  docker run -v <absolute_path_to_data_directory>:/data akshayayadav/overcl-detection-correction run_analysis_subset-F2.sh -c <n>
  ```

