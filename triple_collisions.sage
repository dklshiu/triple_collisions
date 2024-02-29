from sage.modules.misc import gram_schmidt

d=100
Sphere=SphericalDistribution(d)

# Generate a random triple of d-dimensional unit vectors whose pairwise dot products are alpha, beta, gamma 
def rand_triple(alpha,beta,gamma):
   #define the vectors
   a = vector(Sphere.get_random_element())
   b = vector(Sphere.get_random_element())
   c = vector(Sphere.get_random_element())
   #create a list of vectors
   M = [a,b,c]
   # perform Gram Schmidt orthogonalization, inputing matrix M
   # outputs G(an array of vectors) and mu (?)
   G,mu=gram_schmidt(M)
   #with a special structure:
   beta = alpha*gamma + sqrt(1-alpha^2)*(beta-alpha*gamma)/(sqrt(1-alpha^2))
   lam = sqrt(gamma^2 + ((beta-alpha*gamma)^2)/(1-alpha^2))
   bnew = alpha*G[0] + sqrt((1.0-alpha^2))*G[1]
   cnew = gamma*G[0] + (beta - alpha*gamma)*G[1]/sqrt(1.0-alpha^2)
   cnew = cnew + sqrt(1.0-lam^2)*G[2]
   return a,bnew,cnew

# Define a hash function to take the absolute dot product with a fixed set of random vectors and return the m indices with the smallest absolute value (highest values are complementary)
def hash_l(rvecs, target, m):
   dots = { abs(target.dot_product(rvecs[i])):i for i in   range(len(rvecs)) } 
#absolute value of the dot products of the target and each vector in rvecs
   smallest = sorted(dots)[:m] 
#sort in ascending order + keep the first m shortest vectors
   return {dots[z] for z in smallest}

# Run tests for hashing with h random vectors against triple with dot products summing to dot sum, dividing range of dot products into steps
def run_trials(num_trials, h, dot_sum, steps):
   delt = -1.0/steps 
   dot_trips = []

   # Make list of alpha, beta, gamma triples that we wish to test
   num_trips = 0
   for a in range(steps):
      alpha = delt * a
      for b in range(steps):
         beta = delt * b
         gamma = dot_sum - alpha- beta
         if(gamma<0.0 and gamma>-1.0):
            dot_trips += [(alpha,beta,gamma)]
            num_trips += 1

   # For each alpha, beta, gamma make a random triple with that configuration
   oblate_targets = []
   for trip in dot_trips:
      alpha,beta,gamma = trip
      a,b,c = rand_triple(alpha, beta, gamma)
      oblate_targets +=[(a,b,c)]
   # Keep count of collisions for each of our triples
   hist_obl = [0 for j in range(num_trips)] 

   for trials in range(num_trials):
      # Take h fresh random vectors
      # rps = random points on the sphere
      rps = [ vector(Sphere.get_random_element()) for i in range(h)]
      # Hash our oblate triples, counting three-way collisions
      # Hashing converts data into a fixed-length string of letters and numbers
      for k in range(num_trips):
         hash0 = hash_l(rps,oblate_targets[k][0],1)
         hash1 = hash_l(rps,oblate_targets[k][1],1)
         hash2 = hash_l(rps,oblate_targets[k][2],1)
   # Check for collision
         if hash0==hash1==hash2:
            hist_obl[k]+=1

   # Create a list of empirical probabilities using histrogram data
   probs=[]
   for k in range(num_trips):
      print(f"{hist_obl[k]/num_trials*1.0}", dot_trips[k], hist_obl[k])
      probs.append(hist_obl[k]/num_trials*1.0)
   return probs

#Example usage

num_trials = 5000
h = 4

dot_sum = -1.2
steps = 12
p_results = run_trials(num_trials,h,dot_sum,steps)
