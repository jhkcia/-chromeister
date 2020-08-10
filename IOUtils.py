def storeRepresentationMatrix(representationMatrix, location):
    print("Writing matrix")
    with open(location, 'w') as f:
        for i in representationMatrix:
            for j in i:
                f.write("%s " % j)
            f.write("\n")


def loadRepresentationMatrix(location):
    matrix = []
    with open(location, "r") as f:
        for line in f.readlines():
            row = []
            for element in line.split():
                row.append(int(element))
            matrix.append(row)
    return matrix


def storeXYMatrix(representationMatrix, location):
    with open(location, 'w') as f:
        f.write('X Y Count\n')
        for x in range(len(representationMatrix)):
            for y in range(len(representationMatrix[0])):
                if representationMatrix[x][y] > 0:
                    f.write("{} {} {} \n".format(
                        x, y, representationMatrix[x][y]))


def storeEvents(events, location):
    with open(location, 'w') as f:
        for e in events:
            for item in e:
                f.write("%s " % item)
            f.write("\n")
