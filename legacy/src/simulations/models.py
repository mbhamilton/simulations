from django.db import models

# Create your models here.
class Generations(models.Model):
    generation = models.IntegerField()
    
    def __unicode__(self):
        return self.generation