import numpy
from orbkit import grid

from .tools import meshgrid2

def view_with_mayavi(x,y,z,data,is_vectorfield=False,geo_spec=None,datalabels=None,
                     iso_min=1e-4,iso_val=0.01,iso_max=10.0,**kwargs):
  ''' Creates an interactive mayavi dialog showing isosurface plots of the input
  data. 
  
  Components adapted from:
  http://stackoverflow.com/a/1830192
  http://docs.enthought.com/mayavi/mayavi/auto/example_mlab_interactive_dialog.html
  
  **Parameters:**
  
  x,y,z : numpy.ndarray, 1-dim
    Contains the grid.
  data : numpy.ndarray, shape=(len(x), len(y),len(z)) or shape=(N, len(x), len(y),len(z))
    Contains the output data.
  geo_spec : 
    See :ref:`Central Variables` for details.
    If not None, the atom positions will be drawn additionally.
  datalabels : None or list of str
    Contains information about the plotted data with len(datalabels) == len(data).
  '''
  try:
    from enthought.traits.api import HasTraits, Range, Instance, on_trait_change, Bool, Str,List,Button
    from enthought.traitsui.api import View, Item, Group, ListStrEditor,HSplit,VSplit
  except ImportError:
    from traits.api import HasTraits, Range, Instance, on_trait_change, Bool, Str,List,Button
    from traitsui.api import View, Item, Group, ListStrEditor,HSplit,VSplit
  
  try:
    from enthought.mayavi import mlab
    from enthought.mayavi.core.api import PipelineBase
    from enthought.mayavi.core.ui.api import MayaviScene, SceneEditor, MlabSceneModel
  except ImportError:
    from mayavi import mlab
    from mayavi.core.api import PipelineBase
    from mayavi.core.ui.api import MayaviScene, SceneEditor, MlabSceneModel
  
  from copy import deepcopy
  data = numpy.array(data)
  
  if data.ndim == 3:
    data = data[numpy.newaxis]
  elif is_vectorfield and data.ndim == 4:
    data = data[numpy.newaxis]

  if geo_spec is not None and numpy.ndim(geo_spec) == 2:
      geo_spec = [geo_spec for i in range(len(data))]
  
  if datalabels is not None and len(datalabels) < len(data):
    raise ValueError('`datalabels` has to be a list of strings with the same' +
                     'length as `data`.')
  if datalabels is not None:
    datalabels = ['%03d: %s' % (i,j) for i,j in enumerate(datalabels[:len(data)])]
  
  Z,Y,X = meshgrid2(z,y,x)
  
  class MyModel(HasTraits):  
      select  = Range(0, len(data)-1, 0)
      last_select = deepcopy(select)
      iso_value = Range(iso_min, iso_max, iso_val,mode='logslider',
                        label='Iso value' if not is_vectorfield else 'Scale factor')
      if not is_vectorfield:
        opacity = Range(0, 1.0, 1.0)
      else:
        opacity = Range(iso_min, 10*iso_max, 10*iso_val,mode='logslider',
                        label='Range to ')
      show_atoms = Bool(True)
      label = Str()
      available = List(Str)
      available = datalabels
      
      prev_button = Button('Previous')
      next_button = Button('Next')
      
      scene = Instance(MlabSceneModel, ())
      plot_atoms = Instance(PipelineBase)
      plot0 = Instance(PipelineBase)
      
      # When the scene is activated, or when the parameters are changed, we
      # update the plot.
      @on_trait_change('select,iso_value,show_atoms,opacity,label,scene.activated')
      def update_plot(self):        
        if self.plot0 is None:    
          if not is_vectorfield:
            src = mlab.pipeline.scalar_field(X,Y,Z,data[self.select])
            self.plot0 = self.scene.mlab.pipeline.iso_surface(
                        src, contours= [-self.iso_value,self.iso_value], 
                        opacity=self.opacity,colormap='blue-red',
                        vmin=-1e-8,vmax=1e-8)
          else:
            self.plot0 = self.scene.mlab.quiver3d(X,Y,Z,*data[self.select],
                                                  scale_factor=self.iso_value) #flow
          self.plot0.scene.background = (1,1,1)
        elif self.select != self.last_select:
          if not is_vectorfield:
            self.plot0.mlab_source.set(scalars=data[self.select])
          else:
            self.plot0.mlab_source.set(vectors=data[self.select].reshape((3,-1)).T,
                                                  scale_factor=self.iso_value)
          if geo_spec is not None: 
            self.plot_atoms.mlab_source.set(
                                      x=geo_spec[self.select][:,0],
                                      y=geo_spec[self.select][:,1],
                                      z=geo_spec[self.select][:,2])
        if not is_vectorfield:
          self.plot0.contour.contours = [-self.iso_value,self.iso_value]
          self.plot0.actor.property.opacity = self.opacity  
        else:
          #self.plot0.glyph.glyph.range = numpy.array([ 0., self.opacity])
          self.plot0.glyph.glyph.scale_factor = self.iso_value
            
        self.last_select = deepcopy(self.select)
        if datalabels is not None:
          self.label = datalabels[self.select]
        if geo_spec is not None: 
          if self.plot_atoms is None:
            self.plot_atoms = self.scene.mlab.points3d(
                                      geo_spec[self.select][:,0],
                                      geo_spec[self.select][:,1],
                                      geo_spec[self.select][:,2],
                                      scale_factor=0.75,resolution=20)
          self.plot_atoms.visible = self.show_atoms
      
      def _prev_button_fired(self):
        if self.select > 0:
          self.select -= 1
      def _next_button_fired(self):
        if self.select < len(data)-1:
          self.select += 1
      
      
      # The layout of the dialog created
      items = (Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                      height=400, width=600, 
                      show_label=False),)
      items0 = ()
      if len(data) > 1:
        items0 += (Group('select',
                        HSplit(Item('prev_button', show_label=False),
                               Item('next_button', show_label=False)
                              )),)
      items0 += (Group('iso_value', 'opacity','show_atoms'
                        ),)
              
      if datalabels is not None:
        if len(datalabels) > 1:        
            items1 = (Item('available', 
                       editor=ListStrEditor(title='Available Data',editable=False),
                       show_label=False,style='readonly',width=300
                       ),)
            items0 = HSplit(items0,items1)
        items += (Group(Item('label',label='Selected Data',style='readonly', show_label=True),'_'),
               items0,)
      else:
        items += items0
      view = View(VSplit(items[0],items[1:]),
                  resizable=True
                  )
  
  my_model = MyModel()
  my_model.configure_traits()
