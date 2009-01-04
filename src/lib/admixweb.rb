#!/usr/bin/env ruby

$: << File.join(File.expand_path(File.dirname(__FILE__)), %w(.. lib))

require 'admix'

require 'rubygems'
require 'sinatra'

module Admix

  get '/' do
    haml :index
  end

  post '/admix' do |*args|
    begin
      AdmixWrapper.call :loc => params[:loc], :ped => params[:ped]
    rescue AdmixError => e
      throw :halt, [500, "Error: #{e.message}"]
    end
  end

  template :index do <<EOHAML
!!!
%title Admix
%form{ :action => "/admix", :method => "post"}
  %label
    Locus file
    %textarea{ :name => "loc" }
  %label
    Pedigree file
    %textarea{ :name => "ped" }
  %button Submit
EOHAML
  end

  configure do
    set_option :haml, :format => :html5
  end

end
