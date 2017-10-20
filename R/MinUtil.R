SqrtInvCov = function(M)
{
  EigenResult = eigen(as.matrix(M))
  EigenVector = EigenResult$vectors
  EigenValues = abs(EigenResult$values)
  return(EigenVector %*% diag(1/sqrt(EigenValues)) %*% t(EigenVector))
}


Mx = function(M, x)
{
  e = eigen(as.matrix(M))
  return(e$vectors %*% diag(abs(e$values)^x) %*% t(e$vectors))
}

# "%^%" = function(S, pow) with(eigen(S), vectors %*% (abs(values)^pow * t(vectors))) 
# S%^%(-0.5)

mat2ltv = function(mat)
{
  return(mat[upper.tri(mat,diag=TRUE)])
}


ltv2mat = function(vec)
{
  LENGTH = length(vec)
  DIM = round((sqrt(8*LENGTH+1)-1)/2,0)
  if (DIM*(DIM+1)/2 != LENGTH) return(NULL)
  mat = matrix(nrow=DIM, ncol=DIM)
  for (m in 1:DIM) {
    for (n in 1:DIM) {
      k = max(m,n)
      l = min(m,n)
      p = k * (k-1) / 2 + l
      mat[m,n] = vec[p]
    }
  }
  return(mat)
}


mat2utv = function(mat)
{
  return(mat[lower.tri(mat,diag=TRUE)])
}


utv2mat = function(vec)
{
  LENGTH = length(vec)
  DIM = as.integer(round((sqrt(8*LENGTH+1)-1)/2,0))
  if (DIM*(DIM+1)/2 != LENGTH) return(NULL)
  mat = matrix(nrow=DIM, ncol=DIM)
  for (m in 1:DIM) {
    for (n in 1:DIM) {
      k = min(m,n)
      l = max(m,n)
      p = (2*DIM - k + 2)*(k - 1)/2 + l - k + 1
      mat[m,n] = vec[p]
    }
  }
  return(mat)
}


ScaleVar = function(VarMat, dim1)
{
  M1 = chol(VarMat)
  V1 = diag(M1)
  M2 = abs(10 * (M1 - diag(V1, nrow=dim1))) + diag(V1/exp(0.1), nrow=dim1)
  return(t(M2))
}


DesclVar = function(mUCP, mSCL)
{
  nRow = dim(mUCP)[1]
  maT = matrix(nrow=nRow, ncol=nRow)

  for (i in 1:nRow) {
    for (j in 1:nRow) {
      if (i==j) {
        maT[i,j] = exp(mUCP[i,j]) * mSCL[i,j]
      } else if(i > j) {
        maT[i,j] = mUCP[i,j] * mSCL[i,j]
      } else {
        maT[i,j] = 0
      }
    }
  }
  return(maT %*% t(maT))
}


DECN = function(UCP)
{
#Extern: nTheta, nEta, nEps
#        ThetaL[nTheta]==LB, ThetaI[nTheta]==IE, ThetaU[nTheta]==UB
#        alpha, OMscl, SGscl
  uTheta = UCP[1:e$nTheta]
  Theta = exp(uTheta - e$alpha)/(exp(uTheta - e$alpha) + 1)*(e$UB - e$LB) + e$LB
  uvOM  = UCP[(e$nTheta + 1):(e$nTheta + e$nEta*(e$nEta + 1)/2)]
  umOM  = utv2mat(uvOM)
  OM    = DesclVar(umOM, e$OMscl)
  ltvOM = mat2ltv(OM)

  udSG  = UCP[(e$nTheta + e$nEta*(e$nEta + 1)/2 + 1):(e$nTheta + e$nEta*(e$nEta + 1)/2 + e$nEps)]
  umSG  = diag(udSG, nrow=e$nEps)
  SG    = DesclVar(umSG, e$SGscl)
  dgSG  = diag(SG)
  return(c(Theta, ltvOM, dgSG))
}
