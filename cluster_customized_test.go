package cluster

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"testing"

	"github.com/stretchr/testify/assert"
)

type testClusterPoint struct {
	X, Y               float64
	zoom               int
	Id                 int //Index for pint, Id for cluster
	NumPoints          int
	aggregatedClusters []testClusterPoint
}

func TestCluster_Customizer(t *testing.T) {
	points := importData("./testdata/places.json")
	geoPoints := make([]GeoPoint, len(points))
	for i := range points {
		geoPoints[i] = points[i]
	}

	customizer := testCustomizer{}

	c := NewClusterFromCustomizer(customizer)
	c.ClusterPoints(geoPoints, -1)

	northWest := SimplePoint{71.36718750000001, -83.79204408779539}
	southEast := SimplePoint{-71.01562500000001, 83.7539108491127}

	var result []ClusterPoint = c.GetClusters(northWest, southEast, 2)
	for _, p := range result {
		assert.Equal(t, len(p.(testClusterPoint).aggregatedClusters), p.NumberOfPoints()-1, "We get all points except itself")
	}
}

func (cp testClusterPoint) Coordinates() (float64, float64) {
	return cp.X, cp.Y
}

func (cp testClusterPoint) SetCoordinates(X, Y float64) ClusterPoint {
	cp.X = X
	cp.Y = Y
	return cp
}

func (cp testClusterPoint) SetZoom(zoom int) ClusterPoint {
	cp.zoom = zoom
	return cp
}

func (cp testClusterPoint) Zoom() int {
	return cp.zoom
}

func (cp testClusterPoint) SetNumberOfPoints(numPoints int) ClusterPoint {
	cp.NumPoints = numPoints
	return cp
}

func (cp testClusterPoint) NumberOfPoints() int {
	return cp.NumPoints
}

func (cp testClusterPoint) SetID(id int) ClusterPoint {
	cp.Id = id
	return cp
}

func (cp testClusterPoint) ID() int {
	return cp.Id
}

type testCustomizer struct {
}

func (dc testCustomizer) GeoPoint2ClusterPoint(point GeoPoint) ClusterPoint {
	cp := testClusterPoint{aggregatedClusters: make([]testClusterPoint, 0, 10)}
	return cp
}

func (dc testCustomizer) AggregateClusterPoints(point ClusterPoint, aggregated []ClusterPoint, zoom int) ClusterPoint {
	cp := point.(testClusterPoint)
	if cap(cp.aggregatedClusters) < len(cp.aggregatedClusters)+len(aggregated) {
		newAggregated := make([]testClusterPoint, len(cp.aggregatedClusters), max(len(cp.aggregatedClusters)+len(aggregated), 2*cap(cp.aggregatedClusters)))
		copy(cp.aggregatedClusters, newAggregated)
		cp.aggregatedClusters = newAggregated
	}

	for _, p := range aggregated {
		testPoint := p.(testClusterPoint)
		if cap(cp.aggregatedClusters) < len(cp.aggregatedClusters)+1+len(testPoint.aggregatedClusters) {
			newAggregated := make([]testClusterPoint, len(cp.aggregatedClusters), max(len(cp.aggregatedClusters)+1+len(testPoint.aggregatedClusters), 2*cap(cp.aggregatedClusters)))
			cp.aggregatedClusters = newAggregated
		}
		cp.aggregatedClusters = append(cp.aggregatedClusters, testPoint.aggregatedClusters...)
		cp.aggregatedClusters = append(cp.aggregatedClusters, testPoint)
	}
	return cp
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func importData(filename string) []*TestPoint {
	var points = struct {
		Type     string
		Features []*TestPoint
	}{}
	raw, err := ioutil.ReadFile(filename)
	if err != nil {
		fmt.Println(err.Error())
		return nil
	}
	json.Unmarshal(raw, &points)
	//fmt.Printf("Gett data: %+v\n",points)
	return points.Features
}

type TestPoint struct {
	Type       string
	Properties struct {
		//we don't need other data
		Name       string
		PointCount int `json:"point_count"`
	}
	Geometry struct {
		Coordinates []float64
	}
}

func (tp *TestPoint) GetCoordinates() GeoCoordinates {
	return GeoCoordinates{
		Lon: tp.Geometry.Coordinates[0],
		Lat: tp.Geometry.Coordinates[1],
	}
}
