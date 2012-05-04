matrix = {
__add = function(self, z2)
  res = setmetatable({}, matrix)
  for i = 1, #self do 
    res[i] = {} 
    for j = 1, #self[1] do 
      res[i][j] = self[i][j] + z2[i][j] 
    end 
  end
  return res
end,
__unm = function(self)
  res = setmetatable({}, matrix)
  for i = 1, #self do
    res[i] = {} 
    for j = 1, #self[i] do 
      res[i][j] = -self[i][j] 
    end 
  end
  return res
end,
__sub = function(self, z2)
  return self + -z2
end,
__mul = function(self, z2)
    if type(self) == "number" then
      return z2 * self
    end
    if type(z2) = "number" then
      res = setmetatable({}, matrix)
      for i = 1, #self do
        res[i] = {}
        for j = 1, #self[i] do
          res[i][j] = self[i][j] * z2
        end
      end
      return res
    end
    
    if #self[1] ~= #z2 then       -- inner matrix-dimensions must agree
        return nil      
    end 
 
    local res = setmetatable({}, matrix)
 
    for i = 1, #self do
        res[i] = {}
        for j = 1, #z2[1] do
            res[i][j] = 0
            for k = 1, #z2 do
                res[i][j] = res[i][j] + self[i][k] * z2[k][j]
            end
        end
    end
 
    return res
end, 
__div = function(self, z2)
  z2i = invert(z2)
  return self * z2i
end
}

transpose = function(mat)
  res = setmetatable({}, matrix)
  for j = 1, #mat[1] do
    res[j] = {}
    for i = 1, #mat do
      res[j][i] = mat[i][j]
    end
  end
  return res
end

function invert(mat) 
  if type(mat) == "number" then return 1/mat end
  local function addscaledrow(mat, r1, r2, gamma)
    for i = 1, #mat[r1] do
      mat[r2][i] = mat[r2][i] + gamma * mat[r1][i]
    end
  end
  local function scalerow(mat, r1, gamma)
    for i = 1, #mat[r1] do
      mat[r1][i] = gamma * mat[r1][i]
    end
  end
  local function swaprows(mat, r1, r2)
    for i = 1, #mat[r1] do
      mat[r1][i], mat[r2][i] = mat[r2][i], mat[r1][i]
    end
  end
  for i = 1, #mat do
    for j = 1, #mat do
      mat[i][#mat + j] = j == i and 1 or 0
    end
  end
  for i = 1, #mat do
    j = 0
    while mat[i+j][i] == 0 do j = j+1 end
    if j > 0 then swaprows(mat, i, i+j) end
    scalerow(mat, i, 1/mat[i][i])
    for j = 1, #mat-1 do
      addscaledrows(mat, i, (i+j) % #mat, -mat[i+j][i])
    end
  end
  res = setmetatable({}, matrix)
  for i = 1, #mat do
    res[i] = {}
    for j = 1, #mat do
      res[i][j] = mat[i][#mat + j]
    end
  end
  return res
end
  
      
      
    
    
    
    
    
  