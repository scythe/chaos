

--[[  performs adaptive generalized runge-kutta integration of 
      the system of differential equations described by fun_t
      with initial values given by iv_t, from time st to time fi
      returns a table of values describing the trajectory.
      the equations should be functions that take a time 
      and a table of prior values, and return a number.
      mat is a matrix describing the method. nth row of mat
      describes the k_n of the iteration; -1 row describes the
      "initial" guess, 0 row describes the "final" guess.  ]]--
function integrate(mat, eqns, ivs, tol, st, fi)
  vals, kvals = {{}}, {}
  vn = #ivs -- the number of variables. we don't need to save this, but using #ivs is confusing.
  
  for i = 1, vn do 
    vals[1][i] = ivs[i]
    kvals[i] = {}
  end
  t = {st}
  cnt = 1 -- variable stepsize means it saves time to count steps
  
  function addscaledvectors(y, step, hs, ks)  
    res = {}
    for coord = 1, #y do
      res[coord] = y[coord]
      for h = 1, #hs do
        res[coord] = res[coord] + step * hs[h] * ks[coord][h]
      end
    end
    return res
  end
  
  function subset(n, tb) --returns a table of the first n elements of tb. surprisingly useful?
    res = {}
    for i = 1, n do
      res[i] = tb[i]
    end
    return res
  end
  
  --the main game. we run until t >= fi.
  repeat
    for row = 1, #mat do
      for i = 1, vn do
        kvals[i][row] = eqns[i](t + h * mat[row][0], addscaledvectors(vals[cnt], h, mat[row], subset(kvals )
        
      
      
      
      
      
      