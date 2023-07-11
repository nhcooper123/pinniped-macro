# Function to extract data for plotting from BGB
# Modified version of the original BioGeoBEARS results function
# R script that can be found at http://phylo.wikidot.com/biogeobears
# Copyright Nicholas J. Matzke

#-----------------------------------------------
# Get area and state lists
#-----------------------------------------------
extract_results <- function(tree, tipranges, results_object){
  
  tr_pruningwise <- reorder(tree, "pruningwise")
  tips <- 1:length(tr_pruningwise$tip.label)
  nodes <- (length(tr_pruningwise$tip.label) + 1):(length(tr_pruningwise$tip.label) + 
                                                     tr_pruningwise$Nnode)
  
  areas <- getareas_from_tipranges_object(tipranges)
  numareas <- length(areas)
  max_range_size <- 4
  
  numstates <- numstates_from_numareas(numareas = length(areas), 
                                       maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
  
  states_list_areaLetters <- areas_list_to_states_list_new(areas, 
                                                           maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
  
  states_list_0based_index <- rcpp_areas_list_to_states_list(areas, 
                                                             maxareas = max_range_size, include_null_range = results_object$inputs$include_null_range)
  # Set up matrix for node values
  leftright_nodes_matrix <- get_leftright_nodes_matrix_from_results(tr_pruningwise)
  
  # Get marginal probs
  marprobs <- results_object$ML_marginal_prob_each_state_at_branch_bottom_below_node
  
  left_ML_marginals_by_node <- marprobs[leftright_nodes_matrix[, 2], ]
  right_ML_marginals_by_node = marprobs[leftright_nodes_matrix[, 1], ]
  
  # Extract relative probabilities
  relprobs_matrix <- results_object$ML_marginal_prob_each_state_at_branch_top_AT_node
  relprobs_matrix_for_internal_states <- relprobs_matrix[nodes, ]
  
  # Get list of ranges
  ranges_list <- states_list_0based_to_ranges_txt_list(state_indices_0based = states_list_0based_index, 
                                                       areanames = areas)
  statenames <- unlist(ranges_list)
  
  # Extract ML probabilities at nodes
  MLprobs <- get_ML_probs(relprobs_matrix)
  MLstates <- get_ML_states_from_relprobs(relprobs_matrix, statenames, 
                                          returnwhat = "states", if_ties = "takefirst")
  
  # Set up dataframe for plotting
  dd2 <- data.frame(relprobs_matrix_for_internal_states)
  colnames(dd2) <- ranges_list
  dd2$node <- (length(tree$tip.label)+1):((length(tree$edge)/2)+1) 
  dd2$ML <- MLstates[(length(tree$tip.label)+1):((length(tree$edge)/2)+1) ] 
  
  return(dd2)
}