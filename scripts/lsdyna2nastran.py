import re

def writeFloatField(number, length, decimals):
  fmt ='%.' + str(decimals) + 'e'
  # print ('format ' + fmt)
  s = fmt % number
  spaces = ''
  for i in range ((int)(length - len(s))):
    spaces = spaces + ' '
  output = spaces + s
  # print (spaces + s)
  return output


def read_ls_dyna_file(file_path):
    """Reads LS-DYNA keyword file and extracts nodes and triangle elements."""
    nodes = {}
    triangles = []

    with open(file_path, 'r') as f:
        lines = f.readlines()

    in_node_section = False
    in_element_section = False

    for line in lines:
        line = line.strip()
        if line.startswith('$#'):
          continue
        if line.startswith('*NODE'):
            in_node_section = True
            in_element_section = False
            continue
        elif line.startswith('*ELEMENT_SHELL'):
            in_element_section = True
            in_node_section = False
            continue
        elif line.startswith('*'):
            in_node_section = False
            in_element_section = False

        if in_node_section:
            parts = re.split(r'\s+', line)
            if len(parts) >= 4:
                node_id = int(parts[0])
                x, y, z = map(float, parts[1:4])
                nodes[node_id] = (x, y, z)

        if in_element_section:
            parts = re.split(r'\s+', line)
            if len(parts) >= 6:  # LS-DYNA shell elements typically have 4 nodes
                element_id = int(parts[0])
                node_ids = list(map(int, parts[2:6]))
                triangles.append((element_id, node_ids))

    return nodes, triangles


def write_nastran_file(nodes, triangles, file_path):
    """Writes nodes and triangle elements to a NASTRAN small field format file."""
    with open(file_path, 'w') as f:
        f.write("CEND\nBEGIN BULK\n")

        # Write nodes
        for node_id, (x, y, z) in nodes.items():
            f.write(f"{'GRID'.ljust(8)}{str(node_id).rjust(8)}{'0'.rjust(8)}"
                    f"{f'{x:3.4f}'.rjust(8)}{f'{y:3.4f}'.rjust(8)}{f'{z:3.4f}'.rjust(8)}\n")
            #f.write(f"{'GRID'.ljust(8)}{str(node_id).rjust(8)}{'0'.rjust(8)}")
            #f.write(writeFloatField(x,8,4)+writeFloatField(y,8,4)+writeFloatField(z,8,4)+"\n")
        # Write triangles
        for element_id, node_ids in triangles:
            f.write(f"{'CTRIA3'.ljust(8)}{str(element_id).rjust(8)}{'0'.rjust(8)}"
                    f"{str(node_ids[0]).rjust(8)}{str(node_ids[1]).rjust(8)}{str(node_ids[2]).rjust(8)}\n")

        f.write("ENDDATA\n")


# Example usage
ls_dyna_file = 'mesh.k'  # Replace with your LS-DYNA keyword file path
nastran_file = 'mesh.bdf'  # Replace with your output NASTRAN file path

nodes, triangles = read_ls_dyna_file(ls_dyna_file)
write_nastran_file(nodes, triangles, nastran_file)

print(f"Conversion completed. NASTRAN file saved to {nastran_file}.")
