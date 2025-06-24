import utils

file_name = f"TODO_{utils.TAG}.txt"
path = "tmp/"

with open(path + file_name, "a") as file:
  file.write("""
             Figure out and finish the landau fit distributions, include cpu usage times
             Statistics on how to determine if a distribution is better
             check correlations btwn other average DEDx calculation methods\n""")
print(f"File {file_name} created and written to {path}")

