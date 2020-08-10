pro group_min20, inspectra,inback

res=strsplit(inspectra,'.',/extract)

phaname_corr20=res[0]+'_min20.pha'

spawncomando='grppha infile="'+inspectra+'" outfile="'+phaname_corr20+'" chatter=0 comm="chkey BACKFILE "'+inback+'" chkey ANCRFILE /usr/local/swift/caldb41/data/swift/xrt/cpf/arf/swxpc0s6_20010101v013.arf & chkey RESPFILE /usr/local/swift/caldb41/data/swift/xrt/cpf/rmf/swxpc0s6_20010101v013.rmf & systematics 0-1023 0.03 & group min 20 & exit"'



print,spawncomando
spawn,spawncomando

end
