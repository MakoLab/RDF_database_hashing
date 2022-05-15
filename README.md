# Interwoven Hash of Vicious Circle Free Graph

## Requirements

* [networkX](https://networkx.org/)
* [rdflib](https://rdflib.readthedocs.io/en/stable/)

## Usage
```bash
RDF.py [-h] -f -a file
```

where 

* ```-f``` indicates the input serialization format (N-Triples: ```nt```, ```ntriples```, ```n-triples```, Turtle: ```turtle```, ```ttl```, ```n3```, ```notation3```, RDF/XML: ```xml```, ```rdfxml```, and JSON-LD: ```json```, ```json-ld```, ```jsonld```), 
* ```-a``` indicates the hash function (MD5, SHA1, SHA256, SHA512, SHA3 256, SHA 3 512, BLAKE2b, BLAKE2s), and 
* ```-h``` displays help information.

## Alternative implementations

* [Java](https://github.com/MakoLab/ViciousCircleFreeInterwovenHash-Java)
* [C#](https://github.com/MakoLab/ViciousCircleFreeInterwovenHash-CSharp)
