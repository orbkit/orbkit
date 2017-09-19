import orbkit
qc = orbkit.main_read('h2o.molden')

print(qc.mo_spec[[True, True]])
print(qc.mo_spec[[0,1]])
