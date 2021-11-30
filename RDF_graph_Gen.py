###This will be a simple tool for generating somewhat random RDF graphs.


import RDF
import random
import string

Blank_keywords = [
    "Person",
    "Silhouette",
    "Brunette",
    "Creature",
    "Mortal",
    "Human Being",
    "Individual",
    "Character",
    "Personage",
]
Ground_keywords = [
    "Liam",
    "Olivia",
    "Noah",
    "Emma",
    "Oliver",
    "Ava",
    "Elijah",
    "Charlotte",
    "William",
    "Sophia",
    "James",
    "Amelia",
    "Benjamin",
    "Isabella",
    "Lucas",
    "Mia",
    "Henry",
    "Evelyn",
    "Alexander",
    "Harper",
]
Predicate_keywords = [
    "Knows",
    "Enjoys",
    "Likes",
    "Hates",
    "Ignores",
    "Acknowledges",
    "Loves",
]

######################################################
### Inserts a random keyword from a predefined list.
### Distinct lists for blank and ground nodes are
### available.
def insert_random_keyword(blank=False):
    if blank:
        return random.choice(Blank_keywords)
    else:
        return random.choice(Ground_keywords)


######################################################
### Inserts a random predicate from a predefined list.
def insert_random_predicate():
    return random.choice(Predicate_keywords)


######################################################
###
def generate_random_RDF_graph(
    number_of_blanks,
    number_of_ground_nodes,
    number_of_BG_edges,
    number_of_grounded_edges,
    ratio=0.7,
):
    blanks = {}
    ground_nodes = {}
    triplets = []
    for i in range(0, number_of_blanks):
        blanks[i] = RDF.RDF_node("blank" + str(i) + insert_random_keyword(True))
    for i in range(0, number_of_ground_nodes):
        Q = (insert_random_keyword() + "").join(
            random.choices(string.ascii_uppercase + string.digits, k=6)
        ) + str(i)
        ground_nodes[i] = RDF.RDF_node(Q)
    for _ in range(0, number_of_BG_edges):
        a = random.choice(range(0, number_of_blanks - 1))
        b = random.choice(range(a + 1, number_of_blanks))
        triplets.append((blanks[a], insert_random_predicate(), blanks[b]))
    _ = 0
    while _ < number_of_grounded_edges:
        a = random.choice(range(0, number_of_ground_nodes))
        b = ""
        blank = False
        if random.uniform(0, 1) > ratio:
            b = random.choice(range(a, number_of_blanks))
            blank = True
        else:
            b = random.choice(range(0, number_of_ground_nodes))
            while b == a:
                b = random.choice(range(0, number_of_ground_nodes))
        if random.uniform(0, 1) > 0.5:
            if not blank:
                T = (ground_nodes[a], insert_random_predicate(), ground_nodes[b])
            else:
                T = (ground_nodes[a], insert_random_predicate(), blanks[b])
        else:
            if not blank:
                T = (ground_nodes[b], insert_random_predicate(), ground_nodes[a])
            else:
                T = (blanks[b], insert_random_predicate(), ground_nodes[a])
        if T not in triplets:
            _ += 1
            triplets.append(T)
    return RDF.RDF_graph(
        blank_nodes=blanks, real_nodes=ground_nodes, triplets_collection=triplets
    )


####Testing time!####

R = generate_random_RDF_graph(
    number_of_blanks=7,
    number_of_ground_nodes=5,
    number_of_BG_edges=6,
    number_of_grounded_edges=5,
    ratio=0.7,
)
R.to_file(".\\testfiles\\Randomly_generated_graph.txt")
