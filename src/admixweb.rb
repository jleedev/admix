$: << File.join(File.expand_path(File.dirname(__FILE__)), 'lib')

require 'rubygems'
require 'sinatra'

require 'admix'

module Admix

  get '/' do
    haml :index
  end

  post '/admix' do |*args|
    begin
      admix = AdmixWrapper.new :loc => params[:loc], :ped => params[:ped]
      admix.out
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
