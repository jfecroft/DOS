import yaml
with open('tmp.yml', 'r') as f:
        doc = yaml.load(f)
        print doc
