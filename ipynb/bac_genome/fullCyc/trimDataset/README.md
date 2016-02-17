# Dataset for SIPSim validatin

* Day1 of fullCyc HR-SIP dataset
  * pre-fractionated (bulk soil) samples
  * 12C-Con fractionation samples
* Trimming
  * Why
    * To approximate the richness in the reference genome pool being used for simulations
  * Method
    * Just keeping the N-most abundan taxa, where N = number of genomes in reference pool
    * Performing closure operation to assess the the subcompositional abundances 