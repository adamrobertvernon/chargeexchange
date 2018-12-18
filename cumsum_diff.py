import numpy as np

def cumsum_diff(A):
  # Sort A based on first column
  sA = A[np.argsort(A[:,0]),:]

  # Row mask of where each group ends
  row_mask = np.append(np.diff(sA[:,0],axis=0)!=0,[True])

  # Get cummulative summations and then DIFF to get summations for each group
  cumsum_grps = sA.cumsum(0)[row_mask,1:]
  sum_grps = np.diff(cumsum_grps,axis=0)

  # Concatenate the first unique row with its counts
  counts = np.concatenate((cumsum_grps[0,:][None],sum_grps),axis=0)

  # Concatenate the first column of the input array for final output
  out = np.concatenate((sA[row_mask,0][:,None],counts),axis=1)

  return out
