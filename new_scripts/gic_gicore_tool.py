from qgis.core import (
    QgsProcessing,
    QgsProcessingAlgorithm,
    QgsProcessingParameterFeatureSource,
    QgsProcessingParameterString
)

class ExampleAlgorithm(QgsProcessingAlgorithm):

    INPUT_LAYER = 'INPUT_LAYER'
    TEXT_PARAMETER = 'TEXT_PARAMETER'

    def initAlgorithm(self, config=None):
        # Define input parameters
        self.addParameter(
            QgsProcessingParameterFeatureSource(
                self.INPUT_LAYER,
                'Input vector layer',
                [QgsProcessing.TypeVector]
            )
        )
        self.addParameter(
            QgsProcessingParameterString(
                self.TEXT_PARAMETER,
                'Text parameter'
            )
        )

    def processAlgorithm(self, parameters, context, feedback):
        # Retrieve parameter values
        input_layer = self.parameterAsSource(parameters, self.INPUT_LAYER, context)
        text_parameter = self.parameterAsString(parameters, self.TEXT_PARAMETER, context)

        # Perform processing tasks using input_layer and text_parameter
        # ...

        return {}  # Return an empty dictionary as this is a minimal example

    def name(self):
        return 'example_algorithm'

    def displayName(self):
        return 'Example Algorithm'

    def group(self):
        return 'Examples'

    def groupId(self):
        return 'examples'

    def createInstance(self):
        return ExampleAlgorithm()
