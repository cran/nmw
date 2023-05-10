CombDmExPc = function(dm, ex, pc)
{
  ColDm = toupper(colnames(dm))
  ColEx = toupper(colnames(ex))
  ColPc = toupper(colnames(pc))
  colnames(dm) = ColDm
  colnames(ex) = ColEx
  colnames(pc) = ColPc
  
  if (nrow(dm) != length(unique(dm$ID))) stop("The first table should have only one row for each ID!")
  
  ColAll = union(union(ColEx, ColPc), ColDm)
  ToAddEx = setdiff(ColAll, ColEx)
  ToAddPc = setdiff(ColAll, ColPc)

  IDs = sort(unique(pc$ID))
  nID = length(IDs)

  if ("DATE" %in% ColEx) {
    ex$DT = strptime(paste(ex$DATE, ex$TIME), "%Y-%m-%d %H:%M")
    pc$DT = strptime(paste(pc$DATE, pc$TIME), "%Y-%m-%d %H:%M")
  } else {
    ex$DT = ex$TIME
    pc$DT = pc$TIME
  }
  ex = ex[order(ex$ID, ex$DT), ]
  pc = pc[order(pc$ID, pc$DT), ]

  FLAG = rep(T, NROW(ex))
  for (i in 1:nID) {
    cID = IDs[i]
    cDAT = pc[pc$ID == cID, , drop=F]
    cLast = cDAT[NROW(cDAT), "DT"]
    FLAG[ex$ID == cID & ex$DT >= cLast] = F
  }
  ex = ex[FLAG, ]

  ex = cbind(ex, MDV = 1)
  pc = cbind(pc, MDV = 0)

  ToEx = matrix(nrow=NROW(ex), ncol=length(ToAddEx))
  colnames(ToEx) = ToAddEx

  ToPc = matrix(nrow=NROW(pc), ncol=length(ToAddPc))
  colnames(ToPc) = ToAddPc

  Res = rbind(cbind(ex, ToEx), cbind(pc, ToPc)) # Ex first for MDV descending order
  Res = Res[order(Res$ID, Res$DT, Res$MDV), ]
  for (i in 1:nID) {
    cID = IDs[i]
    Res[Res$ID == cID, setdiff(ColDm, "ID")] = dm[dm$ID == cID, setdiff(ColDm, "ID")]
  }
  Res = Res[, c(ColAll, "MDV")]
  rownames(Res) = NULL
  return(Res)
}
