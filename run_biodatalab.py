from assistant.agent import A1

agent = A1(path='./data', llm='gemini-2.5-pro')

agent.go('''Filter sequences with 16-2700 residues in the uniprot_sprot.fasta and save them as file uniprot_sprot_16_2700.fasta''')