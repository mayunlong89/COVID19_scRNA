#date: 2021-04-23
#glb-biotech@zju.edu.cn

##############---1---##################
#1. Hypergeometric test for CXCR6+ memory CD8+T cells

#Severe vs normal

phyper(29,3127,2287,67,lower.tail = T)


#Moderate vs normal

phyper(200,9024,2287,238,lower.tail = F)


#Mild vs normal

phyper(235,8988,2287,273,lower.tail = F)




##############---2---##################
#2. Hypergeometric test for CCR1+CD16+monocytes
#Severe vs normal
phyper(369,3634,1434,416,lower.tail = F)

#Moderate vs normal
phyper(586,6555,1434,633,lower.tail = F)

#Mild vs normal
phyper(162,4081,1434,209,lower.tail = F)


#Druggable genes test

tableR <- matrix(c(40,7,311,136),nrow=2,ncol=2)
chisq.test(tableR)
fisher.test(tableR)


tableR <- matrix(c(22,3,329,140),nrow=2,ncol=2)
chisq.test(tableR)
fisher.test(tableR)

#Hypergeometric test for druggable proteins
phyper(40,351,143,47,lower.tail = F)

#Hypergeometric test for COVID-19 associated druggable proteins
phyper(22,351,143,25,lower.tail = F)


##############---3---##################
#Hypergeometric test for ABO+megakaryocytes
#Severe vs normal
phyper(197,1425,124,203,lower.tail = F)

#Moderate vs normal
phyper(134,1182,124,140,lower.tail = F)

#Mild vs normal
phyper(5,314,124,11,lower.tail = F)


#Hypergeometric test for druggable proteins
phyper(32,424,1558,154,lower.tail = F)

#Hypergeometric test for COVID-19 associated druggable proteins
phyper(26,424,1558,106,lower.tail = F)


#Hypergeometric test for COVID-19 associated druggable proteins
phyper(22,351,143,25,lower.tail = F)



